using System;
using System.IO;
using Sharpmake;

namespace Qaws
{
    // Scalar type fragment for float/double toggle
    [Fragment, Flags]
    public enum ScalarType
    {
        Float  = 1 << 0,  // float (32-bit)
        Double = 1 << 1   // double (64-bit)
    }

    // SIMD fragment for batch evaluation
    [Fragment, Flags]
    public enum SimdMode
    {
        Off  = 1 << 0,
        On   = 1 << 1
    }

    // Custom target with ScalarType and SimdMode
    [Generate]
    public class QawsTarget : Target
    {
        public QawsTarget()
            : base()
        {
        }

        public QawsTarget(Platform platform, DevEnv devEnv, Optimization optimization, ScalarType scalarType, SimdMode simd = SimdMode.Off)
            : base(platform, devEnv, optimization)
        {
            Scalar = scalarType;
            Simd = simd;
        }

        public ScalarType Scalar;
        public SimdMode Simd;
    }

    // Common project base for all Qaws projects
    public abstract class CommonProject : Project
    {
        public string WorkingDir = Path.Combine("[project.SharpmakeCsPath]", "..", "..", "..", "working_dir");

        public CommonProject() : base(typeof(QawsTarget))
        {
            IsFileNameToLower = false;
            IsTargetFileNameToLower = false;
            SourceRootPath = @"[project.SharpmakeCsPath]\..\..\..\src";

            // Add all target permutations
            AddTargets(new QawsTarget(
                Platform.win64,
                DevEnv.vs2022,
                Optimization.Debug | Optimization.Release,
                ScalarType.Float | ScalarType.Double,
                SimdMode.Off | SimdMode.On
            ));
        }

        [Configure()]
        public virtual void ConfigureAll(Configuration conf, QawsTarget target)
        {
            // Name includes scalar type and SIMD: Debug_f32, Release_f64_simd, etc.
            string scalarSuffix = (target.Scalar == ScalarType.Float) ? "f32" : "f64";
            string simdSuffix = (target.Simd == SimdMode.On) ? "_simd" : "";
            conf.Name = "[target.Optimization]_" + scalarSuffix + simdSuffix;

            conf.ProjectFileName = "[project.Name]_[target.DevEnv]_[target.Platform]";
            conf.ProjectPath = Path.Combine("[project.SharpmakeCsPath]", "..", "..", "..", "projects", "[project.Name]");

            conf.IntermediatePath = Path.Combine("[project.SharpmakeCsPath]", "..", "..", "..", "tmp", "[project.Name]", conf.Name);
            conf.TargetPath = WorkingDir;
            // one working_dir for both precisions: f32 binaries get their own names, else an up-to-date f32 build leaves
            // the f64 executable in place
            if (target.Scalar == ScalarType.Float)
                conf.TargetFileSuffix = "_f32";
            conf.TargetLibraryPath = Path.Combine("[project.SharpmakeCsPath]", "..", "..", "..", "tmp", "lib", "win64_" + conf.Name);

            // C11 standard (use compiler option directly)
            conf.AdditionalCompilerOptions.Add("/std:c11");

            // Warning level 4, treat warnings as errors
            conf.AdditionalCompilerOptions.Add("/W4");
            conf.AdditionalCompilerOptions.Add("/WX");

            // Compile as C code (not C++)
            conf.AdditionalCompilerOptions.Add("/TC");

            // Precise floating point for accurate curve calculations
            conf.AdditionalCompilerOptions.Add("/fp:precise");

            // Basic defines
            conf.Defines.Add("NOMINMAX");
            conf.Defines.Add("WIN32");
            conf.Defines.Add("_CRT_SECURE_NO_WARNINGS");

            // Add scalar type define
            if (target.Scalar == ScalarType.Float)
            {
                conf.Defines.Add("QAWS_SCALAR_IS_FLOAT=1");
                // Disable truncation warnings for float builds (double literals to float)
                conf.AdditionalCompilerOptions.Add("/wd4305");  // truncation from 'double' to 'float'
                conf.AdditionalCompilerOptions.Add("/wd4244");  // conversion from 'double' to 'float', possible loss of data
            }
            else
            {
                conf.Defines.Add("QAWS_SCALAR_IS_FLOAT=0");
            }

            // Optimization settings
            if (target.Optimization == Optimization.Debug)
            {
                conf.Defines.Add("_DEBUG");
                conf.Options.Add(Options.Vc.Compiler.Optimization.Disable);
                conf.Options.Add(Options.Vc.Compiler.RuntimeLibrary.MultiThreadedDebug);
            }
            else
            {
                conf.Defines.Add("NDEBUG");
                conf.Options.Add(Options.Vc.Compiler.Optimization.MaximizeSpeed);
                conf.Options.Add(Options.Vc.Compiler.RuntimeLibrary.MultiThreaded);
                conf.Options.Add(Options.Vc.Compiler.Inline.AnySuitable);
            }

            // SIMD batch evaluation support
            if (target.GetFragment<SimdMode>() == SimdMode.On)
            {
                conf.Defines.Add("QAWS_SIMD_AVX2=1");
                conf.AdditionalCompilerOptions.Add("/arch:AVX2");
            }
            else
            {
                // Exclude batch/ sources when SIMD is off
                conf.SourceFilesBuildExcludeRegex.Add(@"[\\/]batch[\\/]");
            }

            // Set working directory for Visual Studio debugging
            conf.VcxprojUserFile = new Configuration.VcxprojUserFileSettings();
            conf.VcxprojUserFile.LocalDebuggerWorkingDirectory = "$(SolutionDir)working_dir";
        }
    }
}
