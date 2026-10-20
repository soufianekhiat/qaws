using System.IO;
using Sharpmake;

namespace Qaws
{
    // Unified test project that runs all tests consecutively
    [Generate]
    public class QawsTestsProject : CommonProject
    {
        public QawsTestsProject()
        {
            Name = "QawsTests";
            SourceRootPath = @"[project.SharpmakeCsPath]\..\..\..\tests\";

            // Include all test files plus the test runner
            SourceFilesExtensions.Add(".c");

            // Exclude standalone debug tools (have their own main())
            SourceFilesExcludeRegex.Add(@"debug_.*\.c$");
        }

        [Configure()]
        public override void ConfigureAll(Configuration conf, QawsTarget target)
        {
            base.ConfigureAll(conf, target);

            // Output type: Console executable
            conf.Output = Configuration.OutputType.Exe;

            // Include paths
            conf.IncludePaths.Add(@"[project.SharpmakeCsPath]\..\..\..\src\qaws");

            // Link against Qaws library
            conf.AddPrivateDependency<QawsLibProject>(target);

            // The suites branch on the precision (sizeof(qaws_scalar), *_F32): constant conditions by design
            conf.AdditionalCompilerOptions.Add("/wd4127");
            // winbase.h under the conforming preprocessor (the SDK's, not ours)
            conf.AdditionalCompilerOptions.Add("/wd5105");        }
    }
}
