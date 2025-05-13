using System.Text.Json;
using System.CommandLine;



public class App
{
    static async Task<int> Main(string[] args)
    {
        var rootCommand = new RootCommand("Calculate distance based on input allele profiles and cluster samples");

        // global options
        var delim = new Option<string>(
            name: "--delimiter",
            description: "character used to separate columns in profile file",
            getDefaultValue: () => "\t"
        )
            .FromAmong(
                "\t",
                ",",
                " "
            );

        var profileFile = new Argument<string>(
            name: "profile",
            description: "profile file"
        );

        var noHeaderOption = new Option<bool>(
            name: "--no-header",
            description: "specify that your profile file has data in the first line",
            getDefaultValue: () => false
        );

        var firstDataFieldOption = new Option<int>(
            name: "--data-field",
            description: "which column in your profile file the allele data begins (e.g., if some columns contain other metadata)",
            getDefaultValue: () => 1
        );

        var sampleFieldOption = new Option<int>(
            name: "--sample-field",
            description: "which column in your profile file contains the sample names",
            getDefaultValue: () => 0
        );

        // upgma options

        var distanceOption = new Option<string>(
            name: "--distance",
            description: "distance measure to use when constructing distance matrix",
            getDefaultValue: () => "absolute"
        )
            .FromAmong(
                "absolute",
                "normalized"
            );

        var maxTreeHeightValue = new Option<float>(
            name: "--max-height",
            description: "Maximum tree height to allow before branch lengths are capped. 0 is unlimited",
            getDefaultValue: () => 0.0F
        );

        var upgma = new Command("upgma", "Cluster samples using UPGMA algorithm"){
            delim,
            noHeaderOption,
            firstDataFieldOption,
            sampleFieldOption,
            distanceOption,
            maxTreeHeightValue,
            profileFile
        };

        rootCommand.AddCommand(upgma);

        upgma.SetHandler((string profile, string delimiter, bool noHeader, int dataStart, int sampleField, string distanceMetric, float treeHeight) =>
            {
                var alleleProfile = ReadProfileFile(profile, delimiter, sampleField, dataStart, noHeader);
                string fullDistance = ParseDistance(distanceMetric);
                DoUPGMA(alleleProfile, fullDistance, treeHeight); 
            },
            profileFile, delim, noHeaderOption,firstDataFieldOption, sampleFieldOption, distanceOption, maxTreeHeightValue
        );

        return await rootCommand.InvokeAsync(args);
    }

    private static string ParseDistance(string distance)
    {
        if (distance == "normalized") {
            return "normalized_allele_differences";
        }
        else {
            return "absolute_allele_differences";
        }
    }

    private static Dictionary<string, List<string>> ReadProfileFile(string filePath, string delim, int keyField, int dataStart, bool noHeader){
        List<string> lines = File.ReadLines(filePath).ToList();
        if (!noHeader) {lines = lines[1..];}
        var profile = new Dictionary<string, List<string>>();
        foreach (var line in lines) {
            string[] fields = line.Split(delim);
            profile.Add(fields[keyField], fields.ToList()[dataStart..]);
        }
        return profile;
    }

    private static void DoUPGMA(Dictionary<string, List<string>> profile, string distanceMetric, float treeHeight){
        var output = DistanceMatrix.FromJson(profile, distance: distanceMetric, max_tree_height: treeHeight, algorithm: "upgma");
        string tree = output.Tree;
        Console.WriteLine(tree);
    }
}
