using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using Sharpmake;

namespace Qaws
{
    // Drives buildsystem/amalgamate.py, which folds the whole tree into a
    // single qaws.h + qaws.c pair for consumers that have no build system.
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
        public static string Python { get; private set; }

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
            Script = Path.Combine(RepoRoot, "buildsystem", "amalgamate.py");
            OutputRoot = Path.Combine(RepoRoot, "tmp", "amalgam");

            Directory.CreateDirectory(Path.Combine(OutputRoot, SimdDir));
            Directory.CreateDirectory(Path.Combine(OutputRoot, NoSimdDir));

            Python = FindPython();
            if (Python == null)
            {
                Builder.Instance.LogWarningLine(
                    "Amalgamation: no Python 3 interpreter found; QawsAmalgam is disabled. "
                    + "Install Python 3 and re-run generate_projects.bat to build it.");
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

        private static string FindPython()
        {
            // "py" is the Windows launcher; the others cover PATH installs.
            string[][] candidates =
            {
                new[] { "py", "-3" },
                new[] { "python3", null },
                new[] { "python", null },
            };

            foreach (string[] candidate in candidates)
            {
                string args = candidate[1] == null ? "--version" : candidate[1] + " --version";
                try
                {
                    if (Execute(candidate[0], args, null) == 0)
                        return candidate[1] == null ? candidate[0] : candidate[0] + " " + candidate[1];
                }
                catch (Exception)
                {
                    // Not on PATH; try the next one.
                }
            }
            return null;
        }

        private static bool Run(string variant, bool simd)
        {
            string arguments = BuildArguments(variant, simd);
            string exe = Python;
            string prefix = "";

            // Python may be "py -3"; split the launcher flag back off.
            int space = exe.IndexOf(' ');
            if (space >= 0)
            {
                prefix = exe.Substring(space + 1) + " ";
                exe = exe.Substring(0, space);
            }

            int code = Execute(exe, prefix + arguments, RepoRoot);
            if (code != 0)
            {
                Builder.Instance.LogWarningLine(
                    "Amalgamation: amalgamate.py failed for the '{0}' variant (exit {1}).", variant, code);
                return false;
            }
            return true;
        }

        private static string BuildArguments(string variant, bool simd)
        {
            string args = string.Format(
                "\"{0}\" --source-dir \"{1}\" --output-dir \"{2}\"",
                Script, SourceDir, Path.Combine(OutputRoot, variant));
            if (simd)
                args += " --simd";
            return args;
        }

        // The command a build step runs to refresh one variant in place.
        public static string RefreshCommand(string variant, bool simd)
        {
            return Python + " " + BuildArguments(variant, simd);
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

            // The generated qaws.c includes "qaws.h" from beside itself.
            conf.IncludePaths.Add(Path.Combine(Amalgamation.OutputRoot, keep));

            // Regenerate before every build so the pair cannot go stale after
            // an edit under src/qaws/.
            conf.EventPreBuild.Add(Amalgamation.RefreshCommand(keep, simd));
            conf.EventPreBuildDescription = "Regenerating the Qaws amalgamation";
        }
    }
}
