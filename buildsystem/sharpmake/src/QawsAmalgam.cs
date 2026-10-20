using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using Sharpmake;

namespace Qaws
{
    // Drives buildsystem/amalgamate.cmake (run with cmake -P), which folds the
    // whole tree into a single qaws.h + qaws.c pair for consumers that have no
    // build system.
    //
    // Two variants are produced, mirroring the SimdMode fragment: the SIMD one
    // also folds in the batch/ sources. Both are generated up front so the
    // files exist by the time Sharpmake scans for sources.
    public static class Amalgamation
    {
        public const string SimdDir = "simd";
        public const string NoSimdDir = "nosimd";

        public static bool Available { get; private set; }
        public static string RepoRoot { get; private set; }
        public static string OutputRoot { get; private set; }
        public static string SourceDir { get; private set; }
        public static string Script { get; private set; }
        public static string CMake { get; private set; }

        static Amalgamation()
        {
            RepoRoot = FindRepoRoot(Directory.GetCurrentDirectory());
            if (RepoRoot == null)
            {
                Builder.Instance.LogWarningLine(
                    "Amalgamation: could not locate the repo root; QawsAmalgam is disabled.");
                return;
            }

            SourceDir = Path.Combine(RepoRoot, "src", "qaws");
            Script = Path.Combine(RepoRoot, "buildsystem", "amalgamate.cmake");
            OutputRoot = Path.Combine(RepoRoot, "tmp", "amalgam");

            Directory.CreateDirectory(Path.Combine(OutputRoot, SimdDir));
            Directory.CreateDirectory(Path.Combine(OutputRoot, NoSimdDir));

            CMake = FindCMake();
            if (CMake == null)
            {
                Builder.Instance.LogWarningLine(
                    "Amalgamation: cmake not found on PATH; QawsAmalgam is disabled. "
                    + "Install CMake and re-run generate_projects.bat to build it.");
                return;
            }

            Available = Run(NoSimdDir, false) && Run(SimdDir, true);
        }

        // Walk up from the working directory until the source tree shows up.
        private static string FindRepoRoot(string start)
        {
            DirectoryInfo dir = new DirectoryInfo(start);
            while (dir != null)
            {
                if (File.Exists(Path.Combine(dir.FullName, "src", "qaws", "qaws.h")))
                    return dir.FullName;
                dir = dir.Parent;
            }
            return null;
        }

        private static string FindCMake()
        {
            try
            {
                if (Execute("cmake", "--version", null) == 0)
                    return "cmake";
            }
            catch (Exception)
            {
                // Not on PATH.
            }
            return null;
        }

        private static bool Run(string variant, bool simd)
        {
            int code = Execute(CMake, BuildArguments(variant, simd), RepoRoot);
            if (code != 0)
            {
                Builder.Instance.LogWarningLine(
                    "Amalgamation: amalgamate.cmake failed for the '{0}' variant (exit {1}).", variant, code);
                return false;
            }
            return true;
        }

        private static string BuildArguments(string variant, bool simd)
        {
            return string.Format(
                "\"-DQAWS_SOURCE_DIR={0}\" \"-DQAWS_OUTPUT_DIR={1}\" -DQAWS_SIMD={2} -P \"{3}\"",
                SourceDir, Path.Combine(OutputRoot, variant), simd ? "ON" : "OFF", Script);
        }

        // The command a build step runs to refresh one variant in place.
        public static string RefreshCommand(string variant, bool simd)
        {
            return CMake + " " + BuildArguments(variant, simd);
        }

        private static int Execute(string exe, string arguments, string workingDir)
        {
            ProcessStartInfo info = new ProcessStartInfo(exe, arguments);
            info.UseShellExecute = false;
            info.RedirectStandardOutput = true;
            info.RedirectStandardError = true;
            info.CreateNoWindow = true;
            if (workingDir != null)
                info.WorkingDirectory = workingDir;

            using (Process process = Process.Start(info))
            {
                process.StandardOutput.ReadToEnd();
                string stderr = process.StandardError.ReadToEnd();
                process.WaitForExit();
                if (process.ExitCode != 0 && !string.IsNullOrEmpty(stderr))
                    Builder.Instance.LogWarningLine("Amalgamation: {0}", stderr.Trim());
                return process.ExitCode;
            }
        }
    }

    // Static library built from the generated qaws.c. It exports the same
    // symbols as Qaws, so QawsTests can be pointed at either one.
    [Generate]
    public class QawsAmalgamProject : CommonProject
    {
        public QawsAmalgamProject()
        {
            Name = "QawsAmalgam";
            SourceRootPath = Amalgamation.OutputRoot;
            SourceFilesExtensions.Add(".c");

            // Both precisions (docs/parameter_convention_plan.md): qaws.c is
            // compiled through two wrappers beside it, one per precision, never
            // directly (as QawsLib's sources).
            if (Amalgamation.OutputRoot == null)
                return;
            foreach (string variant in new string[] { Amalgamation.SimdDir, Amalgamation.NoSimdDir })
            {
                foreach (int bits in new int[] { 32, 64 })
                {
                    string wrapper = Path.Combine(Amalgamation.OutputRoot, variant, "qaws_amalgam_f" + bits + ".c");
                    string text = "/* generated by QawsAmalgam.cs: the amalgamation as f" + bits + " */\n#undef QAWS_SCALAR_IS_FLOAT\n#define QAWS_SCALAR_IS_FLOAT "
                        + (bits == 32 ? "1" : "0") + "\n#include \"qaws.c\"\n";
                    if (!File.Exists(wrapper) || File.ReadAllText(wrapper) != text)
                        File.WriteAllText(wrapper, text);
                    SourceFiles.Add(wrapper);
                }
            }
        }

        [Configure()]
        public override void ConfigureAll(Configuration conf, QawsTarget target)
        {
            base.ConfigureAll(conf, target);

            conf.Output = Configuration.OutputType.Lib;

            bool simd = target.GetFragment<SimdMode>() == SimdMode.On;
            string keep = simd ? Amalgamation.SimdDir : Amalgamation.NoSimdDir;
            string drop = simd ? Amalgamation.NoSimdDir : Amalgamation.SimdDir;

            // Compile only the variant matching this configuration.
            conf.SourceFilesBuildExcludeRegex.Add(@"[\\/]" + drop + @"[\\/]");

            // Both precisions, in every configuration: qaws.c through its two
            // wrappers only; the configuration's scalar only sets what callers'
            // unsuffixed names mean.
            conf.SourceFilesBuildExcludeRegex.Add(@"[\\/]qaws\.c$");
            conf.Defines.Remove("QAWS_SCALAR_IS_FLOAT=1");
            conf.Defines.Remove("QAWS_SCALAR_IS_FLOAT=0");
            conf.Defines.Add("QAWS_SCALAR_IS_FLOAT=0");
            conf.AdditionalCompilerOptions.Add("/wd4305");  // (the f32 objects: truncation from 'double' to 'float')
            conf.AdditionalCompilerOptions.Add("/wd4244");  // (conversion from 'double' to 'float')

            // The generated qaws.c includes "qaws.h" from beside itself.
            conf.IncludePaths.Add(Path.Combine(Amalgamation.OutputRoot, keep));

            // Regenerate before every build so the pair cannot go stale after
            // an edit under src/qaws/.
            conf.EventPreBuild.Add(Amalgamation.RefreshCommand(keep, simd));
            conf.EventPreBuildDescription = "Regenerating the Qaws amalgamation";
        }
    }
}
