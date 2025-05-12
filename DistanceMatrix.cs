using System.Security.Cryptography.X509Certificates;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading.Tasks;

public enum DistMetrics
{
    normalized_allele_differences,
    absolute_allele_differences
}

public enum Algorithms
{
    upgma
}

public class DistanceMatrix
{
    // DistanceMatrix settings
    private static readonly string[] _implementedDistanceMetrics = Enum.GetNames(typeof(DistMetrics));
    private DistMetrics _distanceMetric = DistMetrics.absolute_allele_differences;
    private static readonly string[] _implementedAlgorithms = Enum.GetNames(typeof(Algorithms));
    private Algorithms _algorithm = Algorithms.upgma;

    // Tree settings
    private float _maxTreeHeight = 0;

    // DistanceMatrix data
    private string[] _samples = null;
    private string[][] _profiles = null;
    private float[,] _distmat = null;
    private string _tree = null;
    private string[] _leafOrder = null;



    public string DistanceMetric
    {
        get { return _distanceMetric.ToString(); }
        set 
        {
            if(!Enum.IsDefined(typeof(DistMetrics), value))
            {
                throw new ArgumentException($"Distance metric {value} not implemented. Must be one of {string.Join(", ", _implementedDistanceMetrics)}");
            }
            _distanceMetric = (DistMetrics)Enum.Parse(typeof(DistMetrics), value);
        }
    }

    public string Algorithm
    {
        get { return _algorithm.ToString(); }
        set 
        {
            if(!Enum.IsDefined(typeof(Algorithms), value))
            {
                throw new ArgumentException($"Algorithm {value} not implemented. Must be one of {string.Join(", ", _implementedAlgorithms)}");
            }
            _algorithm = (Algorithms)Enum.Parse(typeof(Algorithms), value);
        }
    }

    public float MaxTreeHeight
    {
        get { return _maxTreeHeight; }
        set 
        {
            if(_distanceMetric == DistMetrics.normalized_allele_differences && value > 1)
            {
                throw new ArgumentException($"MaxTreeHeight of {value} not appropriate for normalized_allele_difference as trees already can't have a height over 1. Should be a value between 0 and 1.");
            }
            if(_distanceMetric == DistMetrics.absolute_allele_differences && value < 1 && value > 0) // between 0 and 1
            {
                throw new ArgumentException($"MaxTreeHeight of {value} not appropriate for absolute_allele_difference as fractional allele differences are not possible. Should be a positive integer value.");
            }
            _maxTreeHeight = value;
        }
    }

    public string[] Samples
    {
        get { return _samples; }
    }

    public float[,] Distmat
    {
        get 
        {
            if(_distmat != null)
            {
                return (float[,])_distmat.Clone();
            }
            SymmetricDistance();
            if(_distmat != null)
            {
                return (float[,])_distmat.Clone();
            }
            throw new InvalidOperationException("Distance matrix could not be calculated.");
        }
    }

    public string Tree
    {
        get
        {
            if(_tree != null)
            {
                return _tree;
            }
            InferTree();
            ReorderDistmat();
            return _tree;
        }
    }

    public string[] LeafOrder
    {
        get
        {
            if(_leafOrder != null)
            {
                return _leafOrder;
            }
            MatchCollection matchList = Regex.Matches(Tree, @"(?<=\(|,)[^(),:]+(?=:)");
            _leafOrder = matchList.Cast<Match>().Select(match => match.Value).ToArray();
            return _leafOrder;
        }
    }

    public void ReadProfileDict(Dictionary<string, List<string>> profileDict)
    {
        var names = new List<string>();
        var profiles = new List<string[]>();
        foreach (var kvp in profileDict)
        {
            if (kvp.Key == "Headers") continue;
            names.Add(kvp.Key);
            profiles.Add(kvp.Value.ToArray());
        }
        _samples = names.Select(n => Regex.Replace(n, "[() ,\":';]", "_")).ToArray();
        _profiles = profiles.ToArray();
    }

    public void SymmetricDistance()
    {
        switch (_distanceMetric)
        {
            case DistMetrics.absolute_allele_differences:
                AbsoluteAlleleDifferences();
                break;
            
            case DistMetrics.normalized_allele_differences:
                NormalizedAlleleDifferences();
                break;
            
            default:
                throw new Exception($"Requested distance metric: {_distanceMetric} is not implemented");
        }
    }

    private void AbsoluteAlleleDifferences()
    {   
        int nSamples = _samples.Length;
        float[,] distmat = new float[nSamples, nSamples];
        
        Parallel.For(0, nSamples, i => {
            for (int j = i+1; j < nSamples; j++)
            {
                string[] a = _profiles[i];
                string[] b = _profiles[j];

                int diffs = 0;
                for (int k = 0; k < a.Length; k++)
                {
                    if (a[k] != b[k] && a[k] != "" && b[k] != "") diffs += 1;
                }

                distmat[i, j] = diffs;
                distmat[j, i] = diffs;
            }
        });
        _distmat = distmat;
    }
    

    private void NormalizedAlleleDifferences()
    {
        int nSamples = _samples.Length;
        float[,] distmat = new float[nSamples, nSamples];
         Parallel.For(0, nSamples, i => {
            for (int j = i+1; j < nSamples; j++)
            {
                string[] a = _profiles[i];
                string[] b = _profiles[j];

                int diffs = 0;
                int compared_loci = 0;
                for (int k = 0; k < a.Length; k++)
                {
                    if (a[k] != "" && b[k] != "") continue;
                    compared_loci += 1;
                    if (a[k] != b[k]) diffs += 1;
                }     
                
                float relative_diff = (float)diffs/compared_loci;

                distmat[i, j] = relative_diff;
                distmat[j, i] = relative_diff;
            }
        });
        _distmat = distmat;
    }

    private void InferTree()
    {
        switch (_algorithm)
        {
            case Algorithms.upgma:
                upgma();
                break;
            
            default:
                throw new Exception($"Requested tree-inference algorithm: {_algorithm} is not implemented");
        }
    }

    private void upgma()
    {
        float[,] distmat = _distmat;
        var nodes = new List<string>(_samples);
        var nodeHeights = new Dictionary<string, float>();

        while (distmat.Length > 1)
        {   
            // Join closest nodes and update tree-building variables
            int[] indices = GetLowestIndex(distmat);
            float[] brlens = GetBanchLengths(distmat, indices, nodes, nodeHeights);
            string node = $"({nodes[indices[0]]}:{brlens[0]:F5},{nodes[indices[1]]}:{brlens[1]:F5})";
            nodes[indices[0]] = node;
            nodes.RemoveAt(indices[1]);
            
            //add node height to dict, capping if appropriate
            float dist = distmat[indices[0], indices[1]];
            if(_maxTreeHeight != 0) // if max height is 0 then the setting is not being used
                {
                    dist = Math.Min(dist, _maxTreeHeight);
                }
            nodeHeights.Add(node, dist/2);

            // Update distmat ready for next round
            float[] avgs = GetAverageDists(distmat, indices);
            distmat = UpdateDistmat(distmat, indices, avgs);

        }
        string tree = $"{nodes[0]};";
        tree = CollapseInternalZeroLengthNodes(tree);
        _tree = tree;
    }

    private int[] GetLowestIndex(float[,] distmat)
    {
        var lowestIdx = new int[2];
        var minValue = 10000.0;
        var dims = distmat.GetLength(0);
        
        for (int i = 0; i < dims; i++)
        {
            for (int j = i+1; j < dims; j++)
            {
                if (distmat[i, j] < minValue)
                {
                    lowestIdx = [i, j];
                    minValue = distmat[i, j];
                }
            }
        }
        return lowestIdx;
    }

    private float[] GetBanchLengths(float[,] distmat, int[] indices, List<string> nodes, Dictionary<string, float> nodeHeights)
    {
        float[] brlens = new float[2];
        int a = indices[0];
        int b = indices[1];
        float dist = distmat[a, b];
        string nodeA = nodes[a];
        string nodeB = nodes[b];
        
        if(_maxTreeHeight != 0) // if max height is 0 then the setting is not being used
        {
            dist = Math.Min(dist, _maxTreeHeight);
        }
        
        float height = dist/2;
        
        if(nodeHeights.ContainsKey(nodeA))
        {
            brlens[0] = height - nodeHeights[nodeA];
        }
        else
        {
            brlens[0] = height;
        }
        if(nodeHeights.ContainsKey(nodeB))
        {
            brlens[1] = height - nodeHeights[nodeB];
        }
        else
        {
            brlens[1] = height;
        }

        return brlens;
    }

    private float[] GetAverageDists(float[,] distmat, int[] idxs)
    /// <summary>
    /// Calculate the mean average in each row of a symmetric matrix for two columns (i and j), but does not compare value [i,j] with [j,i]
    /// </summary>
    /// <param name="distmat">The symmetric matrix to search for values.</param>
    /// <param name="idxs">The two indices to compare</param>
    /// <returns>
    /// Array with length one less that the input array dimensions
    /// </returns>

    {
        var dims = distmat.GetLength(0);
        var a = idxs[0];
        var b = idxs[1];
        var avgs = new List<float>();

        for (int i=0; i<dims; i++)
        {
            if(i == a)
            {
                avgs.Add(0);
            }
            else
            {
                if(i == b)
                {
                    continue;
                }
                else
                {
                    avgs.Add((distmat[a,i] + distmat[b,i])/2);
                }
            }
        }
        return avgs.ToArray();
    }

    private float[,] UpdateDistmat(float[,] distmat, int[] indices, float[] avgs)
    {
        float[,] new_distmat = new float[distmat.GetLength(0)-1, distmat.GetLength(1)-1];
        var to_update = indices[0];
        var to_delete = indices[1];
        var x = 0; // equivalent of i in new_distmat
        for (int i = 0; i < distmat.GetLength(0); i++)
        {
            var y = 0; // equivalent of j in new_distmat
            if(i==to_delete) continue;

            for (int j = 0; j < distmat.GetLength(0); j++)
            {
                if(j==to_delete) continue;

                if(i==to_update)
                {
                    new_distmat[x,y] = avgs[y];
                }
                else
                {
                    if(j==to_update)
                    {
                        new_distmat[x,y] = avgs[x];
                    }
                    else
                    {
                        new_distmat[x,y] = distmat[i,j];
                    }
                }
                y++;
            }
            x++;
        }
        return new_distmat;
    }

    private string CollapseInternalZeroLengthNodes(string old_tree)
    {
        var zero_regex = new Regex(@"\):0\.?0*(?=[,)])");
        var matches = zero_regex.Matches(old_tree);

        // If no matches found then nothing to collapse
        if (matches.Count == 0) return old_tree;

        // Add each closing paren and its associated open paren to a list
        var indices_to_remove = new List<(int Index, int Length)>();

        foreach (Match match in matches)
        {
            int idx = FindOpenParen(old_tree, match.Index);
            indices_to_remove.Add((idx, 1));
            indices_to_remove.Add((match.Index, match.Length));
        }

        indices_to_remove = indices_to_remove.OrderBy(m => m.Index).ToList();

        StringBuilder new_tree = new StringBuilder(indices_to_remove.Count + 1);

        // First, add up to the first paren to remove
        new_tree.Append(old_tree[0..indices_to_remove[0].Index]); 
        for (int i = 0; i < indices_to_remove.Count - 1; i++)
        {
            // Add the string between each set of parens to remove
            int slice_start = indices_to_remove[i].Index + indices_to_remove[i].Length;
            int slice_end = indices_to_remove[(i + 1)].Index;
            new_tree.Append(old_tree[slice_start..slice_end]);
        }

        // Finally, add the remaining string after the last paren to remove
        int final_slice_start = indices_to_remove[^1].Index + indices_to_remove[^1].Length;
        new_tree.Append(old_tree[final_slice_start..^0]);

        return new_tree.ToString();
    }

    private int FindOpenParen(string tree, int closing_paren_idx)
    {
        int paren_count = 0;
        for (int i = closing_paren_idx-1; i>=0; i--)
        {
            char tree_char = tree[i];

            if (tree[i] != '(' && tree[i] != ')') continue; // skip non-paren characters
            
            if (tree[i] == ')') paren_count++; // found a pair of nested parens
            if (tree[i] == '(')
            {
                if(paren_count == 0) return i; // Matching paren found

                paren_count--; // found open paren to nested parens
            }
        }
        // If no open paren was found then something went wrong
        throw new IndexOutOfRangeException($"Ran out of Newick string while trying to find an open parenthesis. This suggests an invalid Newick string has been generated. The closing parenthesis I was trying to pair is at index {closing_paren_idx} in the tree:{tree}");
    }

    private void ReorderDistmat()
    {
        float[,] old_distmat = Distmat;
        float[,] new_distmat = new float[old_distmat.GetLength(0), old_distmat.GetLength(1)];
        List<int> new_indices = new List<int>();
        foreach(string sample in _samples)
        {
            int idx = Array.FindIndex(LeafOrder, x => x == sample);
            new_indices.Add(idx);
        }
        for(int i = 0; i<_samples.Length-1; i++)
        {
            for(int j = i+1; j<_samples.Length; j++)
            {
                int new_i = new_indices[i];
                int new_j = new_indices[j];
                new_distmat[i,j] = old_distmat[new_i, new_j];
                new_distmat[j,i] = old_distmat[new_j, new_i];
            }
        }
        _distmat = new_distmat;
    }

    public static DistanceMatrix FromJson(
        Dictionary<string, List<string>> json,
        string distance = "absolute_allele_differences",
        string algorithm = "upgma",
        float max_tree_height = 200
        )
    {
        DistanceMatrix dmat = new DistanceMatrix
        {
            DistanceMetric = distance,
            Algorithm = algorithm,
            MaxTreeHeight = max_tree_height
        };
        dmat.ReadProfileDict(json);
        dmat.SymmetricDistance();
        return dmat;
    }

    public static DistanceMatrix FromDistanceMatrix(
        string[] samples,
        float[,] distmat,
        string algorithm = "upgma",
        float max_tree_height = 200
        )
    {
        if (distmat.GetLength(0) != distmat.GetLength(1))
        {
            throw new ArgumentException($"Provided distance matrix must have an equal number of rows and columns. The one provided is {distmat.GetLength(0)}X{distmat.GetLength(1)}");
        }
        if (distmat.GetLength(0) != samples.Length)
        {
            throw new ArgumentException($"Provided distance matrix does not have one row and column per sample. You provided a {distmat.GetLength(0)}X{distmat.GetLength(1)} distance matrix and {samples.Length} samples");
        }
        DistanceMatrix dmat = new DistanceMatrix
        {
            _samples = samples,
            _distmat = distmat,
            Algorithm = algorithm,
            MaxTreeHeight = max_tree_height
        };
        return dmat;
    }
}
