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

        var maxTreeHeightValue = new Option<int>(
            name: "--max-height",
            description: "Maximum tree height to allow before branch lengths are capped. 0 is unlimited",
            getDefaultValue: () => 0
        );

        var upgma = new Command("upgma", "Cluster samples using UPGMA algorithm"){
            delim,
            distanceOption,
            maxTreeHeightValue,
            profileFile
        };

        rootCommand.AddCommand(upgma);

        upgma.SetHandler((string profile, string delimiter, string distanceMetric, int treeHeight) =>
            {
                // Console.WriteLine($"profile={profile}, delim={delimiter}, distance={distanceMetric}");
                DoUPGMA(profile, delimiter, distanceMetric, treeHeight); 
            },
            profileFile, delim, distanceOption, maxTreeHeightValue
        );


        return await rootCommand.InvokeAsync(args);

        FileStream json = File.OpenRead(args[0]);
        var profile = JsonSerializer.Deserialize<Dictionary<string, List<string>>>(json);
        var output = DistanceMatrix.FromJson(profile); //, distance: "normalized_allele_differences", max_tree_height: 0);


        // test loading from a distance matrix
        // string[] samples = output.Samples;
        // float[,] dmat = output.Distmat;
        // DistanceMatrix other = DistanceMatrix.FromDistanceMatrix(samples, dmat);

        string tree = output.Tree;
        Console.WriteLine(tree);

        return 0;
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

    private static Dictionary<string, List<string>> ReadProfileFile(string filePath, string delim){
        List<string> lines = File.ReadLines(filePath).ToList();
        var profile = new Dictionary<string, List<string>>();
        foreach (var line in lines) {
            string[] fields = line.Split(delim);
            profile.Add(fields[0], fields.ToList()[1..]);
        }
        return profile;
    }

    private static void DoUPGMA(string profileFile, string delimiter, string distanceOption, int treeHeight){
        string full_distance = ParseDistance(distanceOption);
        Dictionary<string, List<string>> profile = ReadProfileFile(profileFile, delimiter);
        var output = DistanceMatrix.FromJson(profile, distance: full_distance, max_tree_height: treeHeight);
        string tree = output.Tree;
        Console.WriteLine(tree);
    }
}
