# Python-Pybel-Similarity-Search

In this assignment you will use Pybel to perform a similarity search virtual screen. Your script will take two SMILES files as arguments.  The first file will contain the query molecule(s) and the second file will contain the molecules to search.

60% Credit  
Print out the molecules in the search file (the second file argument) that obey Lipinski's rules. Molecules should be printed as canonical SMILES ('can' format).

70% Credit  
Print out each molecule in the search file that matches Lipinski's rules with its  Tanimoto distance to the first molecule in the query file as determined using  Pybel's default fingerprints. The output should display the full SMILES and the Tanimoto coefficient.

80% Credit  
Sort the search molecules by their Tanimoto similarity and print out the top 10 molecules (that match Lipinski) in order of decreasing similarity.

90% Credit  
Compute the average similarity between all the molecules in the query file and the molecules in the search file. Print the top ten Lipinski matching molecules. Note that the query molecules should not be filtered by Lipinski's rules. Make effective use of numpy in performing your calculations.

100% Credit  
Have your script take a similarity threshold as a third argument. Filter out any molecules where the minimum similarity to any query molecule is higher than this threshold. Print the top ten Lipinski matching molecules, sorted by average similarity, after applying this filter.  Along with the molecule, print the minimum, average, and maximum similarity to the query molecules.
