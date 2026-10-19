#ifndef TEST_THREADS_H
#define TEST_THREADS_H

/*
 * Host threads for tests of concurrent evaluation: test_run_threads(n, fn,
 * args) runs fn(args[i]) for i < n, args[0] on the calling thread and the
 * others on threads of their own, and returns when all are done.
 */

#define TEST_MAX_THREADS 16

typedef void (*test_thread_fn)(void* arg);

typedef struct test_thread_arg
{
	test_thread_fn fn;
	void* arg;
} test_thread_arg;

#ifdef _WIN32
#include <windows.h>
#include <process.h>
static unsigned __stdcall test_thread_main(void* p)
{
	test_thread_arg const* a = (test_thread_arg const*)p;
	a->fn(a->arg);
	return 0;
}
#else
#include <pthread.h>
static void* test_thread_main(void* p)
{
	test_thread_arg const* a = (test_thread_arg const*)p;
	a->fn(a->arg);
	return NULL;
}
#endif

static void test_run_threads(unsigned int n, test_thread_fn fn, void* const* args)
{
	test_thread_arg a[TEST_MAX_THREADS];
	unsigned int t;
#ifdef _WIN32
	HANDLE h[TEST_MAX_THREADS];
#else
	pthread_t h[TEST_MAX_THREADS];
#endif
	n = n > TEST_MAX_THREADS ? TEST_MAX_THREADS : n;
	for (t = 0; t < n; t++)
	{
		a[t].fn = fn;
		a[t].arg = args[t];
	}
	for (t = 1; t < n; t++)
#ifdef _WIN32
		h[t] = (HANDLE)_beginthreadex(NULL, 0, test_thread_main, &a[t], 0, NULL);
#else
		pthread_create(&h[t], NULL, test_thread_main, &a[t]);
#endif
	if (n)
		fn(args[0]);
	for (t = 1; t < n; t++)
	{
#ifdef _WIN32
		WaitForSingleObject(h[t], INFINITE);
		CloseHandle(h[t]);
#else
		pthread_join(h[t], NULL);
#endif
	}
}

#endif /* TEST_THREADS_H */
