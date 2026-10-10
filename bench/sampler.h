// A small sampling profiler for the benchmark (Windows only): a thread
// suspends the profiled thread every millisecond, walks its stack, and
// counts self and inclusive samples per function (DbgHelp symbols from the
// PDB). Build with debug information (/Zi, /DEBUG) to get names.

#pragma once

#ifdef _WIN32
#define NOMINMAX
#include <windows.h>
#include <dbghelp.h>
#include <algorithm>
#include <atomic>
#include <cstdio>
#include <map>
#include <set>
#include <string>
#include <thread>
#include <vector>

#pragma comment(lib, "dbghelp.lib")
#pragma comment(lib, "winmm.lib")

class sampler
{
public:
	void start()
	{
		HANDLE proc = GetCurrentProcess();
		DuplicateHandle(proc, GetCurrentThread(), proc, &target_, 0, FALSE, DUPLICATE_SAME_ACCESS);
		SymSetOptions(SYMOPT_DEFERRED_LOADS | SYMOPT_UNDNAME);
		SymInitialize(proc, nullptr, TRUE);
		stop_ = false;
		thread_ = std::thread([this] { run(); });
	}

	void stop_and_report(char const* title, unsigned top = 40)
	{
		stop_ = true;
		thread_.join();
		std::map<std::string, unsigned> self, incl;
		for (unsigned k = 0; k < count_; k++)
		{
			DWORD64 const* st = &frames_[(size_t)k * max_depth];
			std::set<std::string> seen;
			for (unsigned i = 0; i < depth_[k]; i++)
			{
				std::string n = name(st[i]);
				if (i == 0) self[n]++;
				if (seen.insert(n).second) incl[n]++;
			}
		}
		std::vector<std::pair<unsigned, std::string>> vs, vi;
		for (auto& p : self) vs.push_back({ p.second, p.first });
		for (auto& p : incl) vi.push_back({ p.second, p.first });
		std::sort(vs.rbegin(), vs.rend());
		std::sort(vi.rbegin(), vi.rend());
		double n = count_ ? (double)count_ : 1.0;
		std::printf("\n== %s: %zu samples\n-- self\n", title, (size_t)count_);
		for (unsigned i = 0; i < top && i < vs.size(); i++) std::printf("%6.2f%%  %s\n", 100.0 * vs[i].first / n, vs[i].second.c_str());
		std::printf("-- inclusive\n");
		for (unsigned i = 0; i < top && i < vi.size(); i++) std::printf("%6.2f%%  %s\n", 100.0 * vi[i].first / n, vi[i].second.c_str());
		std::fflush(stdout);
		CloseHandle(target_);
		SymCleanup(GetCurrentProcess());
	}

private:
	void run()
	{
		timeBeginPeriod(1);
		while (!stop_ && count_ < max_samples)
		{
			Sleep(1);
			if (SuspendThread(target_) == (DWORD)-1)
				continue;
			CONTEXT ctx;
			memset(&ctx, 0, sizeof(ctx));
			ctx.ContextFlags = CONTEXT_FULL;
			if (GetThreadContext(target_, &ctx))
			{
				/* unwind by the image's unwind tables: no allocation while the thread is stopped */
				DWORD64* st = &frames_[(size_t)count_ * max_depth];
				unsigned d = 0;
				while (d < max_depth && ctx.Rip)
				{
					DWORD64 base = 0;
					PRUNTIME_FUNCTION fe;
					st[d++] = ctx.Rip;
					fe = RtlLookupFunctionEntry(ctx.Rip, &base, nullptr);
					if (!fe)
					{
						/* a leaf: the return address is on top of the stack */
						ctx.Rip = *(DWORD64*)ctx.Rsp;
						ctx.Rsp += 8;
					}
					else
					{
						void* hd = nullptr;
						DWORD64 ef = 0;
						RtlVirtualUnwind(UNW_FLAG_NHANDLER, base, ctx.Rip, fe, &ctx, &hd, &ef, nullptr);
					}
				}
				depth_[count_++] = d;
			}
			ResumeThread(target_);
		}
		timeEndPeriod(1);
	}

	std::string name(DWORD64 a)
	{
		auto it = cache_.find(a);
		if (it != cache_.end()) return it->second;
		char buf[sizeof(SYMBOL_INFO) + 256];
		SYMBOL_INFO* s = (SYMBOL_INFO*)buf;
		memset(buf, 0, sizeof(buf));
		s->SizeOfStruct = sizeof(SYMBOL_INFO);
		s->MaxNameLen = 255;
		DWORD64 disp = 0;
		std::string n = SymFromAddr(GetCurrentProcess(), a, &disp, s) ? std::string(s->Name) : std::string("?");
		cache_[a] = n;
		return n;
	}

	HANDLE target_ = nullptr;
	std::thread thread_;
	std::atomic<bool> stop_{ false };
	static unsigned const max_samples = 20000, max_depth = 48;
	std::vector<DWORD64> frames_ = std::vector<DWORD64>((size_t)max_samples * max_depth);
	std::vector<unsigned> depth_ = std::vector<unsigned>(max_samples);
	unsigned count_ = 0;
	std::map<DWORD64, std::string> cache_;
};
#endif
