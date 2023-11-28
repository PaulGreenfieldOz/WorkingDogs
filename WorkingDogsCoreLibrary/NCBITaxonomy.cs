
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;

namespace WorkingDogsCore
{
    public class NCBITaxonomy
    {
        public static void LoadNamesFile(string namesFN, Dictionary<string, List<int>> nameToTaxId, Dictionary<int, string> taxIdToName)
        {
            StreamReader namesFile = new StreamReader(namesFN);
            char[] vbSeparator = new char[] { '|' };
            bool EOF = false;

            while (!EOF)
            {
                string line = namesFile.ReadLine();
                if (line == null)
                    break;
                if (line == "")
                    continue;

                line = line.Replace("\t", "");
                string[] splitLine = line.Split(vbSeparator);

                int taxId = Convert.ToInt32(splitLine[0]);
                string name = splitLine[1];
                string note = splitLine[2];
                string nameClass = splitLine[3];

                if (nameClass == "scientific name" || nameClass == "synonym")
                {
                    if (!nameToTaxId.ContainsKey(name))
                        nameToTaxId.Add(name, new List<int>());
                    nameToTaxId[name].Add(taxId);

                    if (!taxIdToName.ContainsKey(taxId) && nameClass == "scientific name")
                    {
                        taxIdToName.Add(taxId, name);
                        if (note != "")
                        {
                            nameToTaxId.Add(note, new List<int>());
                            nameToTaxId[note].Add(taxId);
                        }
                    }
                }
            }
            namesFile.Close();

            Console.WriteLine("NCBI: loaded " + nameToTaxId.Count + " names");

        }

        public static void LoadNodeFile(string nodesFN, Dictionary<int, int> taxIdToParent, Dictionary<int, string> taxIdToRank, Dictionary<int, List<int>> taxIdToChildren)
        {
            StreamReader nodesFile = new StreamReader(nodesFN);
            bool EOF = false;
            char[] vbSeparator = new char[] { '|' };

            while (!EOF)
            {
                string line = nodesFile.ReadLine();
                if (line == null)
                    break;
                if (line == "")
                    continue;

                line = line.Replace("\t", "");
                string[] splitLine = line.Split(vbSeparator);

                int taxId = Convert.ToInt32(splitLine[0]);
                int parentTaxId = Convert.ToInt32(splitLine[1]);
                string rank = splitLine[2];

                taxIdToParent.Add(taxId, parentTaxId);
                taxIdToRank.Add(taxId, rank);
                if (!taxIdToChildren.ContainsKey(parentTaxId))
                    taxIdToChildren.Add(parentTaxId, new List<int>());
                taxIdToChildren[parentTaxId].Add(taxId);
            }
            nodesFile.Close();

            Console.WriteLine("NCBI: loaded " + taxIdToParent.Count + " nodes");
        }

        public static void LoadSupplement(string supplementFN, Dictionary<int, int> taxIdToParent, Dictionary<int, string> taxIdToRank, Dictionary<string, List<int>> nameToTaxId, Dictionary<int, string> taxIdToName)
        {
            if (supplementFN == "null")
                return;

            StreamReader supplementFile = new StreamReader(supplementFN);
            int lastUsedTaxId = 100000000;          // decremented on each use
            char[] tabSeparator = new char[] { '\t' };
            bool EOF = false;

            int namesAdded = 0;
            int namesChanged = 0;
            int namesReplaced = 0;
            int namesDeleted = 0;

            while (!EOF)
            {
                string update = supplementFile.ReadLine();
                if (update == null)
                    break;
                if (update == "")
                    continue;
                if (update[0] == '#')
                    continue;

                string[] splitUpdate = update.Split(tabSeparator);

                string nodeName = splitUpdate[0];   // target name
                string rank = splitUpdate[1];       // new rank or delete/rename
                string otherName = splitUpdate[2];  // new name or parent

                if (rank == "delete")
                {
                    if (!nameToTaxId.ContainsKey(nodeName))
                    {
                        Console.WriteLine(nodeName + " not found for delete: " + update);
                        continue;
                    }

                    List<int> taxIds = nameToTaxId[nodeName];
                    if (taxIds.Count > 1)
                    {
                        Console.WriteLine(nodeName + " not deleted. Named mapped to " + taxIds.Count + " taxIds");
                        continue;            
                    }

                    nameToTaxId.Remove(nodeName);

                    if (otherName != "")
                    {
                        List<int> otherTaxIds = nameToTaxId[otherName];
                        if (otherTaxIds.Count > 1 || otherTaxIds[0] != taxIds[0])
                        {
                            Console.WriteLine("delete: taxIds must match. " + nodeName + "(" + taxIds[0] + "), " + otherName + "(" + otherTaxIds[0] + ")");
                            continue;
                        }

                        taxIdToName[taxIds[0]] = otherName;
                    }

                    namesDeleted++;
                    continue;
                }

                if (rank == "rename")
                {
                    if (!nameToTaxId.ContainsKey(nodeName))
                    {
                        Console.WriteLine(nodeName + " not found for rename: " + update);
                        continue;
                    }

                    List<int> taxIds = nameToTaxId[nodeName];
                    if (taxIds.Count != 1)
                    {
                        Console.WriteLine(nodeName + " has " + taxIds.Count + " taxIds. Not renamed: " + update);
                        continue;
                    }

                    if (nameToTaxId.ContainsKey(otherName))
                    {
                        Console.WriteLine(otherName + " already present " + nodeName + " not renamed: " + update);
                        continue;
                    }

                    int taxId = taxIds[0];
                    if (taxIdToName[taxId] == nodeName)
                        taxIdToName[taxId] = otherName;
                    nameToTaxId.Remove(nodeName);
                    nameToTaxId.Add(otherName, taxIds);
                    namesReplaced++;
                    continue;
                }

                if (rank == "synonym")
                {
                    if (nameToTaxId.ContainsKey(nodeName))
                    {
                        Console.WriteLine(nodeName + " already present. " + nodeName + " not added as synonym: " + update);
                        continue;
                    }

                    if (!nameToTaxId.ContainsKey(otherName))
                    {
                        Console.WriteLine(otherName + " not found for rename: " + update);
                        continue;
                    }

                    List<int> taxIds = nameToTaxId[otherName];
                    nameToTaxId.Add(nodeName, taxIds);
                }

                // must be an update or addition

                // check parent name is valid and get its taxId
                if (!nameToTaxId.ContainsKey(otherName))
                {
                    Console.WriteLine(otherName + " (parent) not found for add/change: " + update);
                    continue;
                }
                List<int> parentTaxIds = nameToTaxId[otherName];
                if (parentTaxIds.Count != 1)
                {
                    Console.WriteLine("parent " + otherName + " has " + parentTaxIds.Count + " taxIds. Change/Add not done: " + update);
                    continue;
                }
                int parentTaxId = parentTaxIds[0];

                // does the 'new' name already exist - in which case we're just changing the rank and parent
                // if it does exist, get its taxId, otherwise create new one
                int nodeTaxId = 0;
                bool newNode;

                if (nameToTaxId.ContainsKey(nodeName))
                {
                    List<int> nodeTaxIds = nameToTaxId[nodeName];
                    if (nodeTaxIds.Count != 1)
                    {
                        Console.WriteLine("existing name " + nodeName + " has " + nodeTaxIds.Count + " taxIds. Change not done: " + update);
                        continue;
                    }
                    nodeTaxId = nodeTaxIds[0];
                    newNode = false;
                    namesChanged++;
                }
                else
                {
                    nodeTaxId = lastUsedTaxId;
                    lastUsedTaxId--;
                    newNode = true;
                    namesAdded++;
                }

                // new node so need to add it to the names tables
                if (newNode)
                {
                    nameToTaxId.Add(nodeName, new List<int>());
                    nameToTaxId[nodeName].Add(nodeTaxId);
                    taxIdToName.Add(nodeTaxId, nodeName);

                    // give the new taxName its rank and parent
                    taxIdToParent.Add(nodeTaxId, parentTaxId);
                    taxIdToRank.Add(nodeTaxId, rank);
                }
                else
                {
                    taxIdToParent[nodeTaxId] = parentTaxId;
                    taxIdToRank[nodeTaxId] = rank;
                }
            }

            supplementFile.Close();

            Console.WriteLine("NCBI: added " + namesAdded + ", changed " + namesChanged + ", deleted " + namesDeleted + ", replaced " + namesReplaced);
        }

        public static void SynonymsForNames(Dictionary<string, List<int>> nameToTaxId, Dictionary<int, string> taxIdToName, Dictionary<string, List<string>> nameToSynonyms)
        {
            Dictionary<int, List<string>> namesForTaxid = new Dictionary<int, List<string>>();

            foreach (KeyValuePair<string, List<int>> kvp in nameToTaxId)
            {
                string name = kvp.Key;
                List<int> taxIdsForName = kvp.Value;

                foreach (int taxId in taxIdsForName)
                {
                    if (!namesForTaxid.ContainsKey(taxId))
                        namesForTaxid.Add(taxId, new List<string>());
                    namesForTaxid[taxId].Add(name);
                }
            }

            foreach (KeyValuePair<int, List<string>> kvp in namesForTaxid)
            {
                int taxId = kvp.Key;
                List<string> names = kvp.Value;

                string scientificName = taxIdToName[taxId];

                if (!nameToSynonyms.ContainsKey(scientificName))
                    nameToSynonyms.Add(scientificName, new List<string>());
                foreach (string name in names)
                    nameToSynonyms[scientificName].Add(name);
            }
        }

        public static void NCBITaxonomyFromSpecies(string species, string[] fixedTaxonomySILVA, Dictionary<string, List<int>> nameToTaxId, Dictionary<int, string> taxIdToName, Dictionary<int, string> taxIdToRank,
                                            Dictionary<int, int> taxIdToParent, Dictionary<string, int> rankMappings, HashSet<string> wantedRanks, out List<string[]> fixedTaxonomyNCBINames, out List<int[]> fixedTaxonomyNCBIIds, HashSet<string> reportedClashes)
        {
            int taxId = 0;

            fixedTaxonomyNCBINames = new List<string[]>();
            fixedTaxonomyNCBIIds = new List<int[]>();

            if (species[0] == '\'')
                species = species.Replace("'", "");

            if (species == "unidentified" || species == "uncultured organism" || species.Contains(" symbiont") || species.Contains(" sample"))
                species = "";

            if (species.StartsWith("uncultured "))
                species = species.Substring("uncultured ".Length);

            // does NCBI know about the species name?
            if (species != "" && nameToTaxId.ContainsKey(species))
            {
                taxId = nameToTaxId[species][0];
                if (taxIdToRank[taxIdToParent[taxId]] == "no rank")
                    taxId = 0;
            }

            // try parsing the species name for a genus
            if (taxId == 0)
            {
                int endOfGenus = species.IndexOf(' ');
                if (endOfGenus > 0)
                {
                    string genus = species.Substring(0, endOfGenus);
                    if (genus != "" && nameToTaxId.ContainsKey(genus))
                    {
                        taxId = nameToTaxId[genus][0];
                        if (taxIdToRank[taxId] != "genus")
                            taxId = 0;
                    }
                }
            }

            // didn't find an identifiable species or genus, so try going back up the SILVA taxonomy until we find a NCBI-known name
            if (taxId == 0)
            {
                for (int i = fixedTaxonomySILVA.Length - 1; i >= 0; i--)
                {
                    string taxon = fixedTaxonomySILVA[i];
                    if (taxon != "" && nameToTaxId.ContainsKey(taxon))
                    {
                        taxId = nameToTaxId[taxon][0];
                        if (taxIdToRank[taxId] != "no rank")
                            break;
                    }
                }
            }

            // didn't find any name that NCBI recognised - return empty lists
            if (taxId == 0)
                return;

            // have a name known by NCBI so populate the fixed taxonomy arrays starting from this point
            string taxonName = taxIdToName[taxId];
            foreach (int taxIdInList in nameToTaxId[taxonName])
            {
                string[] fixedTaxonomyNames = new string[rankMappings.Count];
                for (int i = 0; i < rankMappings.Count; i++)
                    fixedTaxonomyNames[i] = "";
                int[] fixedTaxonomyIds = new int[rankMappings.Count];

                taxId = taxIdInList;
                while (taxId != 1)
                {
                    string taxName = taxIdToName[taxId];
                    string rank = taxIdToRank[taxId];

                    if (rankMappings.ContainsKey(rank))
                    {
                        string currentName = fixedTaxonomyNames[rankMappings[rank]];
                        if (currentName != "")
                        {
                            List<int> taxIdsForCurrentName = nameToTaxId[currentName];
                            if (!taxIdsForCurrentName.Contains(taxId))
                            {
                                string clashKey = "NCBI_" + rank + "_" + currentName + "_" + taxName;
                                if (!reportedClashes.Contains(clashKey))
                                {
                                    Console.WriteLine("NCBI: " + rank + " already filled. " + currentName + " -> " + taxName);
                                    reportedClashes.Add(clashKey);
                                }
                            }
                        }
                        fixedTaxonomyNames[rankMappings[rank]] = taxName;
                        fixedTaxonomyIds[rankMappings[rank]] = taxId;
                    }

                    taxId = taxIdToParent[taxId];
                }

                fixedTaxonomyNCBINames.Add(fixedTaxonomyNames);
                fixedTaxonomyNCBIIds.Add(fixedTaxonomyIds);
            }
        }

        public static void TagSILVATaxonomy(string taxonomy, Dictionary<int, string> taxIdToRank, Dictionary<string, List<int>> nameToTaxId, Dictionary<string, int> rankMappings, HashSet<string> wantedRanks, out string[] fixedTaxonomySILVANames, HashSet<string> reportedClashes)
        {
            char[] semiSeparator = new char[] { ';' };
            string[] splitTaxonomy = taxonomy.Split(semiSeparator);
            fixedTaxonomySILVANames = new string[rankMappings.Count];
            for (int i = 0; i < rankMappings.Count; i++)
                fixedTaxonomySILVANames[i] = "";

            // for all the names in the SILVA list... 
            for (int i = 0; i < splitTaxonomy.Length; i++)
            {
                string rank = "";
                string name = splitTaxonomy[i];

                // does NCBI know about this name?
                if (name != "" && nameToTaxId.ContainsKey(name))
                {
                    List<int> taxIdsForName = nameToTaxId[name];

                    // find all the 'wanted' ranks for this name
                    List<string> ranksForName = new List<string>();
                    foreach (int taxId in taxIdsForName)
                    {
                        string rankForTaxId = taxIdToRank[taxId];
                        //if (wantedRanks.Contains(rankForTaxId))
                        ranksForName.Add(rankForTaxId);
                    }

                    rank = ranksForName[0];
                    bool ranksConsistent = true;
                    for (int r = 1; r < ranksForName.Count; r++)
                        if (rank != ranksForName[r])
                            ranksConsistent = false;

                    // ranks are inconsistent so just leave this name out of the 'fixed' set
                    if (!ranksConsistent)
                        continue;

                    // rank both consistent and in the wanted set so record it
                    if (wantedRanks.Contains(rank))
                    {
                        // SILVA:class already filled for Eukaryota;Excavata;Metamonada;Parabasalia;Trichonymphea;Pseudotrichonympha; [1567] Parabasalia -> Trichonymphea
                        string currentName = fixedTaxonomySILVANames[rankMappings[rank]];
                        if (currentName != "")
                        {
                            List<int> taxIdsForCurrentName = nameToTaxId[currentName];
                            if (!taxIdsForName.Contains(taxIdsForCurrentName[0]))
                            {
                                string clashKey = "SILVA_" + rank + "_" + currentName + "_" + name;
                                if (!reportedClashes.Contains(clashKey))
                                {
                                    Console.WriteLine("SILVA: " + rank + " already filled for " + taxonomy + ". " + currentName + " -> " + name);
                                    reportedClashes.Add(clashKey);
                                }
                            }
                        }
                        fixedTaxonomySILVANames[rankMappings[rank]] = name;
                    }
                }
            }
        }
    }
}