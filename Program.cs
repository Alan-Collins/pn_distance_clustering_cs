using System.Text.Json;





public class Test
{
    public static void Main(string[] args)
    {
        FileStream json = File.OpenRead(args[0]);
        var profile = JsonSerializer.Deserialize<Dictionary<string, List<string>>>(json);
        var output = DistanceMatrix.FromJson(profile); //, distance: "normalized_allele_differences", max_tree_height: 0);


        // test loading from a distance matrix
        // string[] samples = output.Samples;
        // float[,] dmat = output.Distmat;
        // DistanceMatrix other = DistanceMatrix.FromDistanceMatrix(samples, dmat);
        
        string tree = output.Tree;
        Console.WriteLine(tree);

    }
}
