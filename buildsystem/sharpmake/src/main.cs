using Sharpmake;

[module: Sharpmake.Include("common.cs")]
[module: Sharpmake.Include("QawsLib.cs")]
[module: Sharpmake.Include("QawsAmalgam.cs")]
[module: Sharpmake.Include("QawsTests.cs")]

namespace Qaws
{
    [Generate]
    public class QawsSolution : Solution
    {
        public QawsSolution() : base(typeof(QawsTarget))
        {
            Name = "Qaws";
            IsFileNameToLower = false;

            // Add targets for both float and double precision
            AddTargets(new QawsTarget(
                Platform.win64,
                DevEnv.vs2022,
                Optimization.Debug | Optimization.Release,
                ScalarType.Float | ScalarType.Double
            ));
        }

        [Configure]
        public void ConfigureAll(Configuration conf, QawsTarget target)
        {
            // Make solution configuration names unique by including scalar type
            string scalarSuffix = (target.Scalar == ScalarType.Float) ? "_f32" : "_f64";
            conf.Name = "[target.Optimization]" + scalarSuffix;

            conf.SolutionFileName = "[solution.Name]_[target.DevEnv]_[target.Platform]";
            conf.SolutionPath = @"[solution.SharpmakeCsPath]\..\..\..";

            // Add library project
            conf.AddProject<QawsLibProject>(target);

            // Add unified test project
            conf.AddProject<QawsTestsProject>(target);

            // The single-file build, when cmake was found on PATH.
            // Building it in every configuration keeps the generated
            // qaws.h/qaws.c pair from silently rotting.
            if (Amalgamation.Available)
                conf.AddProject<QawsAmalgamProject>(target);
        }

        [Main]
        public static void SharpmakeMain(Arguments args)
        {
            args.Generate<QawsSolution>();
        }
    }
}
