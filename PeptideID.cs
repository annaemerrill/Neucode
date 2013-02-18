using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using Coon;
using CoonThermo.IO;
using Coon.Spectra;

namespace OMNE
{
    class PeptideID
    {
        //Constants
        public static Tolerance SILAC = new Tolerance(20.0, ToleranceType.PPM);
        public static Tolerance NEUCODE = new Tolerance(10.0, ToleranceType.PPM);
        public static double INTENSITYCUTOFF = 1.0 / (2.0 * Math.E);
        //public static int NUMCHANNELS = Form1.NUMCHANNELS;
        //public static int NUMISOTOPES = Form1.NUMISOTOPES;

        //Modifications
        public static ChemicalModification CAM = new ChemicalModification("C2H3NO");
        public static ChemicalModification METOX = new ChemicalModification("O");
        
        //Metabolic labels
        public static Mass heavyK = new Mass(8.0142, 8.0142, 8.0142);
        public static Modification heavyLysSILAC = new Modification("Heavy Lys SILAC", heavyK);
        public static double NEUCODELYS = 0.036;
        public static double NEUCODELYS1 = 0.00632;
        public static Mass mediumR = new Mass(6.02013, 6.02013, 6.02013);
        public static Modification mediumArgSILAC = new Modification("Medium Arg SILAC", mediumR);
        public static Mass heavyR = new Mass(10.00827, 10.00827, 10.00827);
        public static Modification heavyArgSILAC = new Modification("Heavy Arg SILAC", heavyR);
        public static Mass averageK = new Mass(8.0322, 8.0322, 8.0322);
        public static Modification averageKNeuCode = new Modification("Average Lys NeuCode", averageK);
        public static Mass heavyKN15 = new Mass(0.997035, 0.997035, 0.997035);
        public static Modification lysPlusOneN15 = new Modification("Heavy Lys+1 NeuCode", heavyKN15);
        public static Mass heavyKC13 = new Mass(1.003355, 1.003355, 1.003355);
        public static Modification lysPlusOneC13 = new Modification("Heavy Lys+1 NeuCode", heavyKC13);
        public static Mass averageK1 = new Mass(1.000195, 1.000195, 1.000195);
        public static Modification averageK1NeuCode = new Modification("Average Lys+1 NeuCode", averageK1);

        //Chemical labels for 6-plex NeuCode
        public static Mass lightmTRAQ = new Mass(140.095, 140.095, 140.095);
        public static Mass lightmTRAQK = new Mass(148.1272, 148.1272, 148.1272);
        public static Modification lightmTRAQNTerm = new Modification("mTRAQ Light N-term", lightmTRAQ);
        public static Modification lightmTRAQLysine = new Modification("mTRAQ Light Lysine", lightmTRAQK);
        public static Modification lightmTRAQTyrosine = new Modification("mTRAQ Light Tyrosine", lightmTRAQ);
        public static Mass mediummTRAQ = new Mass(144.1021, 144.1021, 144.1021);
        public static Mass mediummTRAQK = new Mass(152.1343, 152.1343, 152.1343);
        public static Modification mediummTRAQNTerm = new Modification("mTRAQ Medium N-term", mediummTRAQ);
        public static Modification mediummTRAQLysine = new Modification("mTRAQ Medium Lysine", mediummTRAQK);
        public static Modification mediummTRAQTyrosine = new Modification("mTRAQ Medium Tyrosine", mediummTRAQ);
        public static Mass heavymTRAQ = new Mass(148.1092, 148.1092, 148.1092);
        public static Mass heavymTRAQK = new Mass(156.1414, 156.1414, 156.1414);
        public static Modification heavymTRAQNTerm = new Modification("mTRAQ Heavy N-term", heavymTRAQ);
        public static Modification heavymTRAQLysine = new Modification("mTRAQ Heavy Lysine", heavymTRAQK);
        public static Modification heavymTRAQTyrosine = new Modification("mTRAQ Heavy Tyrosine", heavymTRAQ); 

        //Chemical labels for 12-plex NeuCode
        public static double NEUCODE12 = 0.01264;
        public static Modification MSlight1 = new Modification("light 1", new Mass(431.2210, 431.2210, 431.2210));
        public static Modification MSlight2 = new Modification("light 2", new Mass(431.2336, 431.2336, 431.2336));
        public static Modification MSlight3 = new Modification("light 3", new Mass(431.2462, 431.2462, 431.2462));
        public static Modification MSlight4 = new Modification("light 4", new Mass(431.2588, 431.2588, 431.2588));
        public static Modification MSmedium1 = new Modification("medium 1", new Mass(435.2344, 435.2344, 435.2344));
        public static Modification MSmedium2 = new Modification("medium 2", new Mass(435.2470, 435.2470, 435.2470));
        public static Modification MSmedium3 = new Modification("medium 3", new Mass(435.2596, 435.2596, 435.2596));
        public static Modification MSmedium4 = new Modification("medium 4", new Mass(435.2722, 435.2722, 435.2722));
        public static Modification MSheavy1 = new Modification("heavy 1", new Mass(439.2478, 439.2478, 439.2478));
        public static Modification MSheavy2 = new Modification("heavy 2", new Mass(439.2604, 439.2604, 439.2604));
        public static Modification MSheavy3 = new Modification("heavy 3", new Mass(439.2730, 439.2730, 439.2730));
        public static Modification MSheavy4 = new Modification("heavy 4", new Mass(439.2856, 439.2856, 439.2856));

        //Chemical labels for carbamylation
        public static Modification carbamylLight = new Modification("Carbamyl L", new Mass(44.002849, 44.002849, 44.002849));
        public static Modification carbamylHeavy = new Modification("Carbamyl H", new Mass(44.009168, 44.009168, 44.009168));
        public static Modification carbamylAverage = new Modification("Carbamyl Avg", new Mass(44.0060, 44.0060, 44.0060));

        //Class members
        public Dictionary<ThermoRawFile, List<PeptideSpectralMatch>> PSMs; //Keeps track of PSMs originating from various raw files
        public PeptideSpectralMatch PSM; //Initial PSM
        public int bestPSMScanNumber //Best PSM
        {
            set { }
            get
            {
                int scanNumber = 0;
                double eValue = 1;
                if (PSMs != null && PSMs.Count > 0)
                {
                    foreach (List<PeptideSpectralMatch> psms in PSMs.Values)
                    {
                        if (psms != null && psms[0].EValue < eValue)
                        {
                            eValue = psms[0].EValue;
                            scanNumber = psms[0].ScanNumber;
                        }
                    }
                }
                return scanNumber;
            }
        }
        public LabeledPeak [,] peaks; //Peaks found for the specified number of isotopes for all the quantitative channels
        public double[,] maxIntensity //Finds the peptide's maximum intensity
        {
            get
            {
                double[,] max = null;
                if (allHILACPairs != null && allHILACPairs.Count > 0)
                {
                    max = new double[numChannels, 2];
                    foreach (List<Pair> pairs in allHILACPairs.Values)
                    {
                        foreach (Pair pair in pairs)
                        {
                            for (int i = 0; i < numChannels; i++)
                            {
                                for (int j = 0; j < numIsotopes; j++)
                                {
                                    if (pair.peaks[i,j] != null && pair.peaks[i, j].dNL > max[i, 0])
                                    {
                                        max[i, 0] = pair.peaks[i, j].dNL;
                                        max[i, 1] = pair.rawFile[pair.scanNumber].ScanTime; 
                                    }
                                }
                            }
                        }
                    }
                }
                return max;
            }

        }
        public double[,] maxNoNBCIntensity
        {
            get
            {
                double[,] max = null;
                if (noNBCHILACPairs != null && noNBCHILACPairs.Count > 0)
                {
                    max = new double[numChannels, 1];
                    foreach (List<Pair> pairs in noNBCHILACPairs.Values)
                    {
                        foreach (Pair pair in pairs)
                        {
                            for (int i = 0; i < numChannels; i++)
                            {
                                for (int j = 0; j < numIsotopes; j++)
                                {
                                    if (pair.peaks[i, j] != null && pair.peaks[i, j].dNL > max[i, 0])
                                    {
                                        max[i, 0] = pair.peaks[i, j].dNL;
                                    }
                                }
                            }
                        }
                    }
                }
                return max;
            }
        }
        public double maximumIntensity
        {
            get
            {
                double max = 0;
                if (maxIntensity != null)
                {
                    for (int c = 0; c < numChannels; c++)
                    {
                        if (maxIntensity[c, 0] > max)
                        {
                            max = maxIntensity[c, 0];
                        }
                    }
                }
                return Math.Log10(max);
            }
        }
        public double missingChannelFrequency
        {
            get
            {
                List<Pair> allPairs = new List<Pair>();
                double frequency = -1;
                if (allHILACPairs != null && allHILACPairs.Count > 0)
                {
                    foreach (List<Pair> pairs in allHILACPairs.Values)
                    {
                        foreach (Pair pair in pairs)
                        {
                            allPairs.Add(pair);
                        }
                    }
                    frequency = calculateMissingChannelFrequency(allPairs);
                }
                return frequency;
            }
        }
        public bool coalescenceDetected
        {
            get
            {
                if (coalescedPeakIntensities != null && coalescedPeakIntensities.Count > 0)
                {
                    return true;
                }
                else
                {
                    return false;
                }
            }
        }
        public bool[] quantifiedNBC;
        public bool quantified;
        public List<double> coalescedPeakIntensities;
        public int[] finalQuantified;
        public bool light;
        public bool medium;
        public bool heavy;
        public int numLabels;
        public string sequence;
        public string scanNumbers
        {
            get
            {
                string scanNumbers = "";
                if (PSMs != null && PSMs.Count > 0)
                {
                    foreach (ThermoRawFile rawFile in PSMs.Keys) //For every raw file
                    {
                        List<PeptideSpectralMatch> psms = PSMs[rawFile];
                        foreach (PeptideSpectralMatch PSM in psms) //For every PSM
                        {
                            if (psms.ElementAt(psms.Count - 1) == PSM)
                            {
                                scanNumbers += PSM.ScanNumber;
                            }
                            else
                            {
                                scanNumbers += PSM.ScanNumber + ":";
                            }
                        } //End for every PSM

                        if (PSMs.Keys.ElementAt(PSMs.Keys.Count - 1) != rawFile)
                        {
                            scanNumbers += ";";
                        }
                    } //End for every raw file
                }
                return scanNumbers;
            }
        }
        public string rawFiles
        {
            get
            {
                string rawFiles = "";
                if (PSMs != null && PSMs.Count > 0)
                {
                    for (int i = 0; i < PSMs.Count; i++)
                    {
                        if (i != PSMs.Count() - 1)
                        {
                            rawFiles += PSMs.ElementAt(i).Key + ";";
                        }
                        else
                        {
                            rawFiles += PSMs.ElementAt(i).Key;
                        }
                    }
                }
                return rawFiles;
            }
        }
        public string chargeStates
        {
            get
            {
                string chargeStates = "";
                if (PSMs != null && PSMs.Count > 0)
                {
                    for (int i = 0; i < PSMs.Count; i++)
                    {
                        if (i != PSMs.Count() - 1)
                        {
                            chargeStates += PSMs.ElementAt(i).Value[0].Charge + ";";
                        }
                        else
                        {
                            chargeStates += PSMs.ElementAt(i).Value[0].Charge;
                        }
                    }
                }
                return chargeStates;
            }
        }
        public int numChannels;
        public int numIsotopes;
        public int numClusters;
        public int numIsotopologues;
        public double[,] theoMasses;
        public double[,] adjustedTheoMasses
        {
            get
            {            
                double ppmError = Form1.SYSTEMATICERROR;
                double adjustedMass;
                int channels;
                if (Form1.NEUCODE_SIXPLEX_ARG && Form1.CORRECTPROLINE && prolineCount > 0)
                {
                    channels = numChannels + 2;
                }
                else
                {
                    channels = numChannels;
                }

                double[,] adjustedPrecursorMasses = new double[channels, numIsotopes];
                if (ppmError != 0)
                {
                    for (int i = 0; i < channels; i++)
                    {
                        for (int j = 0; j < numIsotopes; j++)
                        {
                            if (theoMasses[i, 0] > 0)
                            {
                                adjustedMass = ((theoMasses[i, 0] * ppmError) / 1000000.0) + theoMasses[i, 0];
                                adjustedPrecursorMasses[i, j] = adjustedMass + (double)j * Constants.Neutron;
                            }
                            else
                            {
                                adjustedPrecursorMasses[i, j] = 0;
                            }
                        }
                    }
                }
                else
                {
                    for (int i = 0; i < channels; i++)
                    {
                        for (int j = 0; j < numIsotopes; j++)
                        {
                            adjustedMass = theoMasses[i,0];
                            adjustedPrecursorMasses[i, j] = adjustedMass + (double)j * Constants.Neutron;
                        }
                    }
                }
                return adjustedPrecursorMasses;
            }
        }
        public double[,] noNBCTotalIntensity;
        public double[,] totalIntensity;
        public Dictionary<ThermoRawFile, List<Pair>> allHILACPairs;
        public int countAllPairs
        {
            set { }
            get
            {
                int count = 0;
                if (allHILACPairs != null && allHILACPairs.Count > 0)
                {
                    foreach (ThermoRawFile rawFile in allHILACPairs.Keys)
                    {
                        foreach (Pair pair in allHILACPairs[rawFile])
                        {
                            if (pair != null && pair.totalPeakCount > 0)
                            {
                                count++;
                            }
                        }
                    }
                }
                return count;
            }
        }
        public int[] countAllIsotopes
        {
            set { }
            get
            {
                int[] count = null;
                int channelCount;
                int channelIndex;
                if (allHILACPairs != null && allHILACPairs.Count > 0)
                {
                    count = new int[numClusters];
                    foreach (ThermoRawFile rawFile in allHILACPairs.Keys)
                    {
                        List<Pair> pairs = allHILACPairs[rawFile];
                        foreach (Pair pair in pairs)
                        {
                            for (int c = 0; c < numClusters; c++)
                            {
                                for (int j = 0; j < numIsotopes; j++)
                                {
                                    channelCount = 0;
                                    channelIndex = c * numIsotopologues;
                                    for (int i = channelIndex; i < channelIndex + numIsotopologues; i++)
                                    {
                                        if (pair.peaks[i, j] != null)
                                        {
                                            channelCount++;
                                        }
                                    }
                                    if ((Form1.NEUCODE_ARG_PROLINECONVERSION && channelCount > 0) || channelCount == numIsotopologues)
                                    {
                                        count[c]++;
                                    }
                                }
                            }
                        }
                    }
                }
                return count;
            }
        }
        public Dictionary<ThermoRawFile, List<Pair>> noNBCHILACPairs;
        public int countNoNBCPairs
        {
            set { }
            get
            {
                int count = 0;
                if (noNBCHILACPairs != null && noNBCHILACPairs.Count > 0)
                {
                    foreach (ThermoRawFile rawFile in noNBCHILACPairs.Keys)
                    {
                        foreach (Pair pair in noNBCHILACPairs[rawFile])
                        {
                            if (pair != null && pair.totalPeakCount > 0)
                            {
                                count++;
                            }
                        }
                    }
                    return count;
                }
                return count;
            }
        }
        public int[] countNoNBCIsotopes
        {
            set { }
            get
            {
                int[] count = null;
                int channelCount;
                int channelIndex;
                if (noNBCHILACPairs != null && noNBCHILACPairs.Count > 0)
                {
                    count = new int[numClusters];
                    foreach (ThermoRawFile rawFile in noNBCHILACPairs.Keys)
                    {
                        List<Pair> pairs = noNBCHILACPairs[rawFile];
                        foreach (Pair pair in pairs)
                        {
                            for (int c = 0; c < numClusters; c++)
                            {
                                for (int j = 0; j < numIsotopes; j++)
                                {
                                    channelCount = 0;
                                    channelIndex = c * numIsotopologues;
                                    for (int i = channelIndex; i < channelIndex + numIsotopologues; i++)
                                    {
                                        if (pair.peaks[i, j] != null)
                                        {
                                            channelCount++;
                                        }
                                    }
                                    if (channelCount == numIsotopologues)
                                    {
                                        count[c]++;
                                    }
                                }
                            }
                        }
                    }
                }
                return count;
            }
        }
        public Dictionary<ThermoRawFile, PeptideSpectralMatch> bestPSM;
        public Dictionary<ThermoRawFile, int[]> firstLastSN;
        public int lysineCount;
        public int arginineCount;
        public int prolineCount;
        public double[,] heavyToLightRatioSum;
        public Range spacingRange
        {
            get
            {
                Range spacing;
                double theoSpacing;
                if (numLabels > 0)
                {
                    if (Form1.NEUCODE)
                    {
                        if (numIsotopologues == 4)
                        {
                            theoSpacing = (double)numLabels * NEUCODE12;
                            spacing = new Range(theoSpacing, new Tolerance(0.007, ToleranceType.DA));
                        }
                        else if (Form1.NEUCODE_DUPLEX_LYS1 || Form1.NEUCODE_DUPLEX_CARBAMYL)
                        {
                            theoSpacing = (double)numLabels * NEUCODELYS1;
                            spacing = new Range(theoSpacing, new Tolerance(0.002, ToleranceType.DA));
                        }
                        else
                        {
                            theoSpacing = (double)lysineCount * NEUCODELYS;
                            spacing = new Range(theoSpacing, new Tolerance(0.010, ToleranceType.DA));
                        }
                        
                    }
                    else
                    {
                        theoSpacing = (double)lysineCount * heavyK.MonoisotopicMass;
                        spacing = new Range(theoSpacing, new Tolerance(0.020, ToleranceType.DA));
                    }
                    return spacing;
                }
                else
                {
                    return null;
                }
            }
        }
        public Range rTRange;
        public Range fullScanRange;
        public List<int> fullScanList;
        public Tolerance firstSearchRange;        

        //Constructor for duplex HILAC or SILAC
        public PeptideID(int scanNumber, int charge, double theoMass, double eValue, string sequence, ThermoRawFile rawFile, string mods)
        {
            List<PeptideSpectralMatch> psms;
            int[] SN;
            Peptide peptide;
            Peptide checkLight;
            Peptide checkMedium;
            Peptide checkHeavy;
            double checkLightError;
            double checkMediumError;
            double checkHeavyError;
            double theoMassCorrection;
            
            //Initialize all peptide properties
            this.sequence = sequence;
            countLysinesArgininesProlines();
            numChannels = Form1.NUMCHANNELS;
            numIsotopes = Form1.NUMISOTOPES;
            numIsotopologues = Form1.NUMISOTOPOLOGUES;
            numClusters = Form1.NUMCLUSTERS;

            if (Form1.NEUCODE)
            {
                if (numIsotopologues == 4)
                {
                    numLabels = lysineCount + 1;
                    for (int i = 0; i < mods.Length; i++)
                    {
                        if (mods[i].Equals('Y'))
                        {
                            numLabels++;
                        }
                    }
                    //firstSearchRange = new Tolerance((3 * NEUCODE12 * (double)(numLabels) * 0.75), ToleranceType.DA);
                    firstSearchRange = new Tolerance(200.0, ToleranceType.PPM);
                }
                else if (Form1.NEUCODE_DUPLEX_LYS1)
                {
                    numLabels = lysineCount;
                    firstSearchRange = new Tolerance(150.0, ToleranceType.PPM);
                }
                else if (Form1.NEUCODE_DUPLEX_CARBAMYL)
                {
                    numLabels = lysineCount + 1;
                    firstSearchRange = new Tolerance(150.0, ToleranceType.PPM);
                }
                else
                {
                    numLabels = lysineCount;
                    firstSearchRange = new Tolerance((NEUCODELYS * (double)lysineCount * 0.75), ToleranceType.DA);
                }
            }
            else
            {
                if (Form1.NEUCODE_ARG_PROLINECONVERSION)
                {
                    numLabels = arginineCount;
                    firstSearchRange = new Tolerance(50.0, ToleranceType.PPM);
                }
                else
                {
                    numLabels = lysineCount;
                    firstSearchRange = new Tolerance(50.0, ToleranceType.PPM);
                }
            }

            theoMasses = new double[numChannels, 1];

            //Initialize all identification data members
            PSMs = new Dictionary<ThermoRawFile, List<PeptideSpectralMatch>>();
            PSMs.Add(rawFile, new List<PeptideSpectralMatch>());
            PSM = new PeptideSpectralMatch(scanNumber, charge, sequence, null, eValue);
            PSMs.TryGetValue(rawFile, out psms);
            psms.Add(PSM);
            bestPSM = new Dictionary<ThermoRawFile, PeptideSpectralMatch>();
            bestPSM.Add(rawFile, PSM);
            firstLastSN = new Dictionary<ThermoRawFile, int[]>();
            SN = new int[2];
            firstLastSN.Add(rawFile, SN);
            firstLastSN.TryGetValue(rawFile, out SN);
            SN[0] = PSM.ScanNumber;
            SN[1] = PSM.ScanNumber;
            
            //Initialize all quantitation data members
            countNoNBCPairs = 0;
            countNoNBCIsotopes = new int[numClusters];
            countAllPairs = 0;
            countAllIsotopes = new int[numClusters];
            noNBCTotalIntensity = new double[numChannels, numIsotopes + 1];
            totalIntensity = new double[numChannels, numIsotopes + 1];
            if (Form1.NEUCODE)
            {
                heavyToLightRatioSum = new double[(numIsotopologues - 1) * numClusters, 1];
            }
            else
            {
                heavyToLightRatioSum = new double[numChannels / 2, 1];
            }
            allHILACPairs = new Dictionary<ThermoRawFile, List<Pair>>();
            allHILACPairs.Add(rawFile, new List<Pair>());
            noNBCHILACPairs = new Dictionary<ThermoRawFile, List<Pair>>();
            noNBCHILACPairs.Add(rawFile, new List<Pair>()); 
            
            //Set theoretical masses of each channel
            //When multiple isotopic clusters (i.e., separated by more than a Da) are present, verify which cluster the PSM belongs to before setting theo mass
            if (numLabels > 0)
            {                
                peptide = new Peptide(sequence, mods);
                //peptide.SetFixedModifications('R', new Modification("empty R", new Mass(0.0, 0.0, 0.0)));
                //peptide.SetFixedModifications('C', new Modification("empty C", new Mass(0.0, 0.0, 0.0)));
                //peptide.SetFixedModifications('K', new Modification("empty K", new Mass(0.0, 0.0, 0.0)));
                peptide.SetFixedModifications('C', CAM);
             
                if (Form1.NEUCODE)
                {
                    if (Form1.NEUCODE_DUPLEX_LYS1)
                    {
                        theoMassCorrection = (NEUCODELYS1 / 2.0) * (double)lysineCount;
                        checkLight = new Peptide(peptide);
                        checkLight.SetFixedModifications('K', averageK1NeuCode);
                        theoMass = checkLight.Mass.MonoisotopicMass;
                        theoMasses[0, 0] = theoMass - theoMassCorrection;
                        theoMasses[1, 0] = theoMass + theoMassCorrection;
                    }

                    if (Form1.NEUCODE_DUPLEX_CARBAMYL)
                    {
                        theoMassCorrection = (NEUCODELYS1 / 2.0) * (double)numLabels;
                        checkLight = new Peptide(peptide);
                        if (sequence[0].Equals('K'))
                        {
                            bool nonKFound = false;
                            int counter = 1;

                            while (!nonKFound && counter < sequence.Length)
                            {
                                if (!sequence[counter].Equals('K'))
                                {
                                    nonKFound = true;
                                }
                                counter++;
                            }

                            checkLight.SetFixedModification(counter + 1, carbamylAverage);
                        }
                        else
                        {
                            checkLight.SetFixedModification(1, carbamylAverage);
                        }
                        checkLight.SetFixedModifications('K', carbamylAverage);
                        theoMass = checkLight.Mass.MonoisotopicMass;
                        theoMasses[0, 0] = theoMass - theoMassCorrection;
                        theoMasses[1, 0] = theoMass + theoMassCorrection;
                    }

                    //NeuCode SILAC duplex
                    if (Form1.NEUCODE_DUPLEX_LYS)
                    {
                        theoMassCorrection = (NEUCODELYS / 2.0) * (double)lysineCount;
                        checkLight = new Peptide(peptide);
                        checkLight.SetFixedModifications('K', averageKNeuCode);
                        //checkLight.SetFixedModifications('R', heavyArgSILAC);
                        theoMass = checkLight.Mass.MonoisotopicMass;
                        theoMasses[0, 0] = theoMass - theoMassCorrection;
                        theoMasses[1, 0] = theoMass + theoMassCorrection;
                    }

                    if (Form1.NEUCODE_4PLEX_LIGHT)
                    {
                        checkLight = new Peptide(peptide);
                        //double mass = checkLight.Mass.MonoisotopicMass;
                        if (sequence[0].Equals('K'))
                        {
                            bool nonKFound = false;
                            int counter = 1;

                            while (!nonKFound && counter < sequence.Length)
                            {
                                if (!sequence[counter].Equals('K'))
                                {
                                    nonKFound = true;
                                }
                                counter++;
                            }

                            checkLight.SetFixedModification(counter + 1, MSlight1);
                        }
                        else
                        {
                            checkLight.SetFixedModification(1, MSlight1);
                        }
                        checkLight.SetFixedModifications('K', MSlight1);

                        theoMasses[0, 0] = checkLight.Mass.MonoisotopicMass;
                        theoMasses[1, 0] = theoMasses[0, 0] + (numLabels * NEUCODE12);
                        theoMasses[2, 0] = theoMasses[1, 0] + (numLabels * NEUCODE12);
                        theoMasses[3, 0] = theoMasses[2, 0] + (numLabels * NEUCODE12);
                    }
                    
                    if (Form1.NEUCODE_4PLEX_HEAVY)
                    {
                        checkHeavy = new Peptide(peptide);
                        if (sequence[0].Equals('K'))
                        {
                            bool nonKFound = false;
                            int counter = 1;

                            while (!nonKFound && counter < sequence.Length)
                            {
                                if (!sequence[counter].Equals('K'))
                                {
                                    nonKFound = true;
                                }
                                counter++;
                            }

                            checkHeavy.SetFixedModification(counter + 1, MSheavy1);
                        }
                        else
                        {
                            checkHeavy.SetFixedModification(1, MSheavy1);
                        }
                        checkHeavy.SetFixedModifications('K', MSheavy1);
                        
                        theoMasses[0, 0] = checkHeavy.Mass.MonoisotopicMass;
                        theoMasses[1, 0] = theoMasses[0, 0] + (numLabels * NEUCODE12);
                        theoMasses[2, 0] = theoMasses[1, 0] + (numLabels * NEUCODE12);
                        theoMasses[3, 0] = theoMasses[2, 0] + (numLabels * NEUCODE12);
                    }

                    //NeuCode 12-plex
                    if (Form1.NEUCODE_12PLEX)
                    {                        
                        checkLight = new Peptide(peptide);
                        checkMedium = new Peptide(peptide);
                        checkHeavy = new Peptide(peptide);

                        if (sequence[0].Equals('K'))
                        {
                            bool nonKFound = false;
                            int counter = 1;

                            while (!nonKFound && counter < sequence.Length)
                            {
                                if (!sequence[counter].Equals('K'))
                                {
                                    nonKFound = true;
                                }
                                counter++;
                            }

                            checkLight.SetFixedModification(counter + 1, MSlight1);
                            checkMedium.SetFixedModification(counter + 1, MSmedium1);
                            checkHeavy.SetFixedModification(counter + 1, MSheavy1);
                        }
                        else
                        {
                            checkLight.SetFixedModification(1, MSlight1);
                            checkMedium.SetFixedModification(1, MSmedium1);
                            checkHeavy.SetFixedModification(1, MSheavy1);
                        }
                        checkLight.SetFixedModifications('K', MSlight1);
                        checkMedium.SetFixedModifications('K', MSmedium1);
                        checkHeavy.SetFixedModifications('K', MSheavy1);

                        //Set theo masses
                        theoMasses[0, 0] = checkLight.Mass.MonoisotopicMass;
                        theoMasses[1, 0] = theoMasses[0,0] + (numLabels * NEUCODE12);
                        theoMasses[2, 0] = theoMasses[1, 0] + (numLabels * NEUCODE12);
                        theoMasses[3, 0] = theoMasses[2, 0] + (numLabels * NEUCODE12);
                        theoMasses[4, 0] = checkMedium.Mass.MonoisotopicMass;
                        theoMasses[5, 0] = theoMasses[4, 0] + (numLabels * NEUCODE12);
                        theoMasses[6, 0] = theoMasses[5, 0] + (numLabels * NEUCODE12);
                        theoMasses[7, 0] = theoMasses[6, 0] + (numLabels * NEUCODE12);
                        theoMasses[8, 0] = checkHeavy.Mass.MonoisotopicMass;
                        theoMasses[9, 0] = theoMasses[8, 0] + (numLabels * NEUCODE12);
                        theoMasses[10, 0] = theoMasses[9, 0] + (numLabels * NEUCODE12);
                        theoMasses[11, 0] = theoMasses[10, 0] + (numLabels * NEUCODE12);
                    }

                    //NeuCode SILAC w/ mTRAQ
                    if (Form1.NEUCODE_SIXPLEX_MTRAQ)
                    {
                        theoMassCorrection = (NEUCODELYS / 2.0) * (double)lysineCount;
                        
                        checkLight = new Peptide(peptide);
                        checkLight.SetFixedModifications('K', lightmTRAQLysine);
                        checkLight.SetFixedModification(1, lightmTRAQNTerm);

                        checkMedium = new Peptide(peptide);
                        checkMedium.SetFixedModifications('K', mediummTRAQLysine);
                        checkMedium.SetFixedModification(1, mediummTRAQNTerm);

                        checkHeavy = new Peptide(peptide);
                        checkHeavy.SetFixedModifications('K', heavymTRAQLysine);
                        checkHeavy.SetFixedModification(1, heavymTRAQNTerm);

                        checkLightError = Tolerance.GetError(theoMass, checkLight.Mass.MonoisotopicMass, ToleranceType.PPM);
                        checkMediumError = Tolerance.GetError(theoMass, checkMedium.Mass.MonoisotopicMass, ToleranceType.PPM);
                        checkHeavyError = Tolerance.GetError(theoMass, checkHeavy.Mass.MonoisotopicMass, ToleranceType.PPM);

                        //Set as mTRAQ version with smallest mass error
                        if (Math.Abs(checkLightError) <= Math.Abs(checkMediumError) && Math.Abs(checkLightError) <= Math.Abs(checkHeavyError)) //Identification is "light"
                        {
                            //Initialize and update "light"
                            light = true;
                            //mTRAQLight = new PeptideID(scanNumber, charge, checkLight.Mass.MonoisotopicMass, eValue, rawFile, this, checkMedium.Mass.MonoisotopicMass, checkHeavy.Mass.MonoisotopicMass);
                        }
                        else if (Math.Abs(checkMediumError) <= Math.Abs(checkHeavyError) && Math.Abs(checkMediumError) <= Math.Abs(checkLightError)) //Identification is "medium"
                        {
                            //Initialize and update "medium"
                            medium = true;
                            //mTRAQMedium = new PeptideID(scanNumber, charge, checkMedium.Mass.MonoisotopicMass, eValue, rawFile, this, checkLight.Mass.MonoisotopicMass, checkHeavy.Mass.MonoisotopicMass);
                        }
                        else //Identification is "heavy"
                        {
                            //Initialize and update "heavy"
                            heavy = true;
                            //mTRAQHeavy = new PeptideID(scanNumber, charge, checkHeavy.Mass.MonoisotopicMass, eValue, rawFile, this, checkLight.Mass.MonoisotopicMass, checkMedium.Mass.MonoisotopicMass);
                        }

                        //Set theo masses
                        theoMasses[0, 0] = checkLight.Mass.MonoisotopicMass - theoMassCorrection;
                        theoMasses[1, 0] = checkLight.Mass.MonoisotopicMass + theoMassCorrection;
                        theoMasses[2, 0] = checkMedium.Mass.MonoisotopicMass - theoMassCorrection;
                        theoMasses[3, 0] = checkMedium.Mass.MonoisotopicMass + theoMassCorrection;
                        theoMasses[4, 0] = checkHeavy.Mass.MonoisotopicMass - theoMassCorrection;
                        theoMasses[5, 0] = checkHeavy.Mass.MonoisotopicMass + theoMassCorrection;
                    }
                    
                    if (Form1.NEUCODE_SIXPLEX_ARG)
                    {
                        theoMassCorrection = (NEUCODELYS / 2.0) * (double)lysineCount;
                        peptide.SetFixedModifications('K', averageKNeuCode);
                        
                        if (arginineCount == 0)
                        {
                            theoMass = peptide.Mass.MonoisotopicMass;
                            theoMasses[0, 0] = theoMass - theoMassCorrection;
                            theoMasses[1, 0] = theoMass + theoMassCorrection;
                        }
                        else
                        {
                            checkLight = new Peptide(peptide);

                            checkMedium = new Peptide(peptide);
                            checkMedium.SetFixedModifications('R', mediumArgSILAC);

                            checkHeavy = new Peptide(peptide);
                            checkHeavy.SetFixedModifications('R', heavyArgSILAC);

                            checkLightError = Tolerance.GetError(theoMass, checkLight.Mass.MonoisotopicMass, ToleranceType.PPM);
                            checkMediumError = Tolerance.GetError(theoMass, checkMedium.Mass.MonoisotopicMass, ToleranceType.PPM);
                            checkHeavyError = Tolerance.GetError(theoMass, checkHeavy.Mass.MonoisotopicMass, ToleranceType.PPM);

                            //Set as mTRAQ version with smallest mass error
                            if (Math.Abs(checkLightError) <= Math.Abs(checkMediumError) && Math.Abs(checkLightError) <= Math.Abs(checkHeavyError)) //Identification is "light"
                            {
                                //Initialize and update "light"
                                light = true;
                                //mTRAQLight = new PeptideID(scanNumber, charge, checkLight.Mass.MonoisotopicMass, eValue, rawFile, this, checkMedium.Mass.MonoisotopicMass, checkHeavy.Mass.MonoisotopicMass);
                            }
                            else if (Math.Abs(checkMediumError) <= Math.Abs(checkHeavyError) && Math.Abs(checkMediumError) <= Math.Abs(checkLightError)) //Identification is "medium"
                            {
                                //Initialize and update "medium"
                                medium = true;
                                //mTRAQMedium = new PeptideID(scanNumber, charge, checkMedium.Mass.MonoisotopicMass, eValue, rawFile, this, checkLight.Mass.MonoisotopicMass, checkHeavy.Mass.MonoisotopicMass);
                            }
                            else //Identification is "heavy"
                            {
                                //Initialize and update "heavy"
                                heavy = true;
                                //mTRAQHeavy = new PeptideID(scanNumber, charge, checkHeavy.Mass.MonoisotopicMass, eValue, rawFile, this, checkLight.Mass.MonoisotopicMass, checkMedium.Mass.MonoisotopicMass);
                            }        

                            //Set theo masses
                            theoMasses[0, 0] = checkLight.Mass.MonoisotopicMass - theoMassCorrection;
                            theoMasses[1, 0] = checkLight.Mass.MonoisotopicMass + theoMassCorrection;
                            theoMasses[2, 0] = checkMedium.Mass.MonoisotopicMass - theoMassCorrection;
                            theoMasses[3, 0] = checkMedium.Mass.MonoisotopicMass + theoMassCorrection;
                            theoMasses[4, 0] = checkHeavy.Mass.MonoisotopicMass - theoMassCorrection;
                            theoMasses[5, 0] = checkHeavy.Mass.MonoisotopicMass + theoMassCorrection;

                            //Add 1 proline to heavy mass of proline-containing peptides
                            if (Form1.CORRECTPROLINE && prolineCount > 0)
                            {
                                checkHeavy.SetFixedModification(1, new Modification("heavy proline", new Mass(6.01381, 6.01381, 6.01381)));
                                theoMasses[6, 0] = checkHeavy.Mass.MonoisotopicMass - theoMassCorrection;
                                theoMasses[7, 0] = checkHeavy.Mass.MonoisotopicMass + theoMassCorrection;
                            }
                        }
                    }
                }

                if (Form1.NEUCODE_ARG_PROLINECONVERSION)
                {
                    checkLight = new Peptide(peptide);
                    checkLight.SetFixedModifications('R', mediumArgSILAC);

                    //checkMedium = new Peptide(checkLight); // 1 proline 
                    //checkHeavy = new Peptide(checkLight); // 2 proline

                    //checkMedium.SetFixedModifications('P', new Modification("heavy proline", new Mass(6.01381, 6.01381, 6.01381)));
                    //checkMedium.SetFixedModification(1, new Modification("heavy proline", new Mass(15.0504, 15.0504, 15.0504)));
                    //checkHeavy.SetFixedModification(1, new Modification("heavy proline", new Mass(10.0336, 10.0336, 10.0336)));

                    theoMasses[0, 0] = checkLight.Mass.MonoisotopicMass;
                    //theoMasses[1, 0] = checkMedium.Mass.MonoisotopicMass;
                }

                //Traditional SILAC
                if (Form1.SILAC_DUPLEX_LYS)
                {
                    checkLight = new Peptide(peptide);
                    checkHeavy = new Peptide(checkLight);
                    checkHeavy.SetFixedModifications('K', heavyLysSILAC);
                    checkLightError = Tolerance.GetError(theoMass, checkLight.Mass.MonoisotopicMass, ToleranceType.PPM);
                    checkHeavyError = Tolerance.GetError(theoMass, checkHeavy.Mass.MonoisotopicMass, ToleranceType.PPM);

                    if (Math.Abs(checkLightError) <= Math.Abs(checkHeavyError))
                    {
                        light = true;
                    }
                    else
                    {
                        heavy = true;
                    }

                    //Set theo masses
                    theoMasses[0, 0] = checkLight.Mass.MonoisotopicMass;
                    theoMasses[1, 0] = checkHeavy.Mass.MonoisotopicMass;
                }
            }
        }

        //Constructor for 6-plex HILAC child (identified)
        /*public PeptideID(int scanNumber, int charge, double theoMass, double eValue, ThermoRawFile rawFile, PeptideID parent, double theoMass1, double theoMass2)
        {
            this.parent = parent;
            this.lysineCount = parent.lysineCount;
            this.sequence = parent.sequence;

            //Update identified
            averageTheoMass = theoMass;

            if (parent.light) //Light mTRAQ identification
            {
                //Initialize "medium" & "heavy"
                parent.mTRAQMedium = new PeptideID(charge, theoMass1, rawFile, parent);
                parent.mTRAQHeavy = new PeptideID(charge, theoMass2, rawFile, parent);           
            }
            else if (parent.medium) //Medium mTRAQ identification
            {
                //Initialize "light" & "heavy"
                parent.mTRAQLight = new PeptideID(charge, theoMass1, rawFile, parent);
                parent.mTRAQHeavy = new PeptideID(charge, theoMass2, rawFile, parent);
            }
            else //Heavy mTRAQ identification
            {
                //Initialize "light" & "medium"
                parent.mTRAQLight = new PeptideID(charge, theoMass1, rawFile, parent);
                parent.mTRAQMedium = new PeptideID(charge, theoMass2, rawFile, parent);
            }
        }

        //Constructor for 6-plex HILAC child (predicted)
        public PeptideID(int charge, double theoMass, ThermoRawFile rawFile, PeptideID parent)
        {
            this.parent = parent;
            this.sequence = parent.sequence;
            this.lysineCount = parent.lysineCount; 
            averageTheoMass = theoMass;
        }*/

        //Constructor for 6-plex HILAC parent
        /*public PeptideID(int scanNumber, int charge, double theoMass, double eValue, string sequence, ThermoRawFile rawFile, string mods)
        {
            //Initialize all parent members
            //orderedScanNumbers = new Dictionary<ThermoRawFile, SortedList<int, double>>();
            //orderedScanNumbers.Add(rawFile, new SortedList<int, double>());
            //orderedScanNumbers[rawFile].Add(scanNumber, eValue);
            PSM = new PeptideSpectralMatch(scanNumber, charge, sequence, null, eValue);
            bestPSM = new Dictionary<ThermoRawFile, PeptideSpectralMatch>();
            bestPSM.Add(rawFile, PSM);
            firstLastSN = new Dictionary<ThermoRawFile, int[]>();
            firstLastSN.Add(rawFile, new int[2]);
            firstLastSN[rawFile][0] = PSM.ScanNumber;
            firstLastSN[rawFile][1] = PSM.ScanNumber;
            PSMs = new Dictionary<ThermoRawFile, List<PeptideSpectralMatch>>();
            PSMs.Add(rawFile, new List<PeptideSpectralMatch>());
            PSMs[rawFile].Add(PSM);
            this.sequence = sequence;
            countLysines();
            light = false;
            medium = false;
            heavy = false;
            maxIntensity = new double[Form1.NUMCHANNELS, 1];
            firstSearchRange = new Tolerance(0.375 * ((NEUCODELYS * (double)lysineCount) / (double)charge), ToleranceType.DA);

            //Figure out which version was identified
            Peptide peptide = new Peptide(sequence, mods); 
            peptide.SetFixedModifications('C', CAM);

            Peptide checkLight = new Peptide(peptide);
            Peptide checkMedium = new Peptide(peptide);
            Peptide checkHeavy = new Peptide(peptide);

            //Try setting identified peptide as "light", "medium", or "heavy" & check error from theoretical mass
            checkLight.SetFixedModifications('K', lightmTRAQLysine);
            checkLight.SetFixedModification(1, lightmTRAQNTerm);
            
            checkMedium.SetFixedModifications('K', mediummTRAQLysine);
            checkMedium.SetFixedModification(1, mediummTRAQNTerm);

            checkHeavy.SetFixedModifications('K', heavymTRAQLysine);
            checkHeavy.SetFixedModification(1, heavymTRAQNTerm);

            double checkLightError = Tolerance.GetError(theoMass, checkLight.Mass.MonoisotopicMass, ToleranceType.PPM);
            double checkMediumError = Tolerance.GetError(theoMass, checkMedium.Mass.MonoisotopicMass, ToleranceType.PPM);
            double checkHeavyError = Tolerance.GetError(theoMass, checkHeavy.Mass.MonoisotopicMass, ToleranceType.PPM);

            //Set as mTRAQ version with smallest mass error
            if (Math.Abs(checkLightError) <= Math.Abs(checkMediumError) && Math.Abs(checkLightError) <= Math.Abs(checkHeavyError)) //Identification is "light"
            {
                //Initialize and update "light"
                light = true;
                mTRAQLight = new PeptideID(scanNumber, charge, checkLight.Mass.MonoisotopicMass, eValue, rawFile, this, checkMedium.Mass.MonoisotopicMass, checkHeavy.Mass.MonoisotopicMass);      
            }
            else if (Math.Abs(checkMediumError) <= Math.Abs(checkHeavyError) && Math.Abs(checkMediumError) <= Math.Abs(checkLightError)) //Identification is "medium"
            {
                //Initialize and update "medium"
                medium = true;
                mTRAQMedium = new PeptideID(scanNumber, charge, checkMedium.Mass.MonoisotopicMass, eValue, rawFile, this, checkLight.Mass.MonoisotopicMass, checkHeavy.Mass.MonoisotopicMass);
            }
            else //Identification is "heavy"
            {
                //Initialize and update "heavy"
                heavy = true;
                mTRAQHeavy = new PeptideID(scanNumber, charge, checkHeavy.Mass.MonoisotopicMass, eValue, rawFile, this, checkLight.Mass.MonoisotopicMass, checkMedium.Mass.MonoisotopicMass);
            }        
        }*/

        //Counts the number of lysines in the peptide
        public void countLysinesArgininesProlines()
        {
            lysineCount = 0;
            arginineCount = 0;
            prolineCount = 0;
            string peptid = sequence;
            char[] residues = peptid.ToCharArray();
            for (int i = 0; i < residues.Length; i++)
            {
                if (residues[i] == 'K' || residues[i] == 'k')
                {
                    lysineCount++;
                }
                else if (residues[i] == 'R' || residues[i] == 'r')
                {
                    arginineCount++;
                }
                else if (residues[i] == 'P' || residues[i] == 'p')
                {
                    prolineCount++;
                }
            }
        }

        public LabeledPeak largestPeak (double theoMass, ThermoRawFileScan currentScan, Tolerance range, ThermoRawFile rawFile)
        {
            PeptideSpectralMatch best;
            LabeledPeak largest;
            //LabeledPeak closest;
            Range mZRange;

            int charge;
            double dNLIntensity;
            double SN = Form1.MINIMUMSN;

            best = PSMs[rawFile].ElementAt(0);
            charge = best.Charge;

            double theoMZ = Mass.MZFromMass(theoMass, charge);

            mZRange = new Range(theoMZ, range);

            if (theoMass == 0)
            {
                return null;
            }
            else
            {                
                largest = (LabeledPeak)currentScan.Spectrum.GetLargestPeak(mZRange);
                //closest = (LabeledPeak)currentScan.Spectrum.GetNearestPeak(theoMZ);
                //double error = Math.Abs(Tolerance.GetError(closest.MZ, theoMZ, ToleranceType.PPM));
                //Check closest peak
                /*if (error <= 5.0)
                {
                    if (closest.Charge > 0 && closest.Charge != charge)
                    {
                        closest = null;
                    }
                    else if (closest.SN < Form1.MINIMUMSN)
                    {
                        closest = null;
                    }
                }
                else
                {
                    closest = null;
                }
                return closest;
                /*if (largest == null && closest != null)
                {
                    if (closest.SN >= SN)
                    {
                        dNLIntensity = closest.Intensity * rawFile.GetInjectionTime(currentScan.ScanNum); //Calculate dNL intensity
                        closest.Intensity = dNLIntensity;
                        return closest;
                    }
                }
                else if (closest != null && largest.MZ != closest.MZ)
                {
                    dNLIntensity = closest.Intensity * rawFile.GetInjectionTime(currentScan.ScanNum); //Calculate dNL intensity
                    closest.Intensity = dNLIntensity;
                    return closest;
                }*/
                if (largest != null)
                {
                    if (largest.Charge > 0 && largest.Charge != charge) //Ignore peaks that are the wrong charge state while considering those with unknown charge
                    {
                        largest = null;
                    }

                    else if (largest.SN < SN) //Ignore low S/N peaks
                    {
                        largest = null;
                    }

                    else
                    {
                        //dNLIntensity = largest.Intensity * rawFile.GetInjectionTime(currentScan.ScanNum); //Calculate dNL intensity
                        //largest.dNL = dNLIntensity;
                        return largest;
                    }
                }
            }
            return null;
        }

        public void largestPeakPPM (double[] theoMasses, ThermoRawFileScan currentScan, Tolerance range, ThermoRawFile rawFile, List<PrecursorPPM> ppms)
        {
            PeptideSpectralMatch best;
            Range mZRange;

            int charge;
            double SN = Form1.MINIMUMSN;

            best = PSMs[rawFile].ElementAt(0);
            charge = best.Charge;
            double theoMass = 0;

            for (int i = 0; i < theoMasses.Length; i++)
            {
                theoMass += theoMasses[i];
            }

            theoMass = theoMass / theoMasses.Length;

            double theoMZ = Mass.MZFromMass(theoMass, charge);
            mZRange = new Range(theoMZ, range);

            if (theoMass == 0)
            {
                return;
            }
            else
            {
                List<IPeak> peaks = null;
                List<LabeledPeak> sortedPeaks = null;
                List<LabeledPeak> topPeaks = null;
                if (currentScan.Spectrum.TryGetPeaks(out peaks, mZRange))
                {
                    sortedPeaks = new List<LabeledPeak>();
                    foreach (IPeak peak in peaks)
                    {
                        LabeledPeak newPeak = (LabeledPeak)peak;
                        if (newPeak.SN < 5)
                        {
                            newPeak = null;
                        }
                        if (newPeak != null && newPeak.Charge != 0 && newPeak.Charge != charge)
                        {
                            newPeak = null;
                        }
                        if (newPeak != null)
                        {
                            sortedPeaks.Add(newPeak);
                        }
                    }
                    sortedPeaks.Sort(LabeledPeak.sortIntensityDescending());

                    if (numIsotopologues == 4)
                    {
                        if (sortedPeaks.Count > 3)
                        {
                            topPeaks = new List<LabeledPeak>();
                            topPeaks.Add(sortedPeaks[0]);
                            topPeaks.Add(sortedPeaks[1]);
                            topPeaks.Add(sortedPeaks[2]);
                            topPeaks.Add(sortedPeaks[3]);
                            topPeaks.Sort(LabeledPeak.sortMZAscending());

                            PrecursorPPM ppm;

                            for (int i = 0; i < theoMasses.Length; i++)
                            {
                                ppm = new PrecursorPPM(charge, this.sequence, best.EValue, Tolerance.GetError(topPeaks[i].MZ, Mass.MZFromMass(theoMasses[i], charge), ToleranceType.PPM));
                                ppms.Add(ppm);
                            }
                        }
                    }
                    if (numIsotopologues == 2)
                    {
                        if (sortedPeaks.Count > 1)
                        {
                            topPeaks = new List<LabeledPeak>();
                            topPeaks.Add(sortedPeaks[0]);
                            topPeaks.Add(sortedPeaks[1]);
                            topPeaks.Sort(LabeledPeak.sortMZAscending());

                            PrecursorPPM ppm;

                            for (int i = 0; i < theoMasses.Length; i++)
                            {
                                ppm = new PrecursorPPM(charge, this.sequence, best.EValue, Tolerance.GetError(topPeaks[i].MZ, Mass.MZFromMass(theoMasses[i], charge), ToleranceType.PPM));
                                ppms.Add(ppm);
                            }
                        }
                    }
                }
                return;
            }
        }

        /*public void extractNoNBCPairs(ThermoRawFile rawFile)
        {
            List<Pair> pairs;
            Pair pair;
            int numChannels = Form1.NUMCHANNELS;
            int numIsotopes = Form1.NUMISOTOPES;

            bool scanFound = false;

            allHILACPairs.TryGetValue(rawFile, out pairs);
            for (int i = 0; i < pairs.Count(); i++)
            {
                scanFound = false;
                pair = pairs[i];

                for (int j = 0; j < numChannels; j += 2)
                {
                    for (int k = 0; k < numIsotopes; k++)
                    {
                        if (pair.peaks[j, k] != null && pair.peaks[j + 1, k] != null)
                        {
                            countNoNBCIsotopes++;
                            lightNoNBCTotalIntensity[(j/2), k] += pair.peaks[j, k].Intensity;
                            lightNoNBCTotalIntensity[(j/2), numIsotopes] += pair.peaks[j, k].Intensity;
                            heavyNoNBCTotalIntensity[(j/2), k] += pair.peaks[j+1, k].Intensity;
                            heavyNoNBCTotalIntensity[(j/2), numIsotopes] += pair.peaks[j+1, k].Intensity;
                            scanFound = true;
                        }
                    }
                }
                if (scanFound)
                {
                    countNoNBCPairs++;
                }
            }
        }*/

        /* Calculates the MS1 scan range to use for quantification
         * Includes all PSMs, extending above and below according to the entered retention time window
         */
        public void calculateScanRange(ThermoRawFile rawFile, double rTWindow)
        {
            PeptideSpectralMatch best = PSMs[rawFile].ElementAt(0);
            List<PeptideSpectralMatch> allPSMs = PSMs[rawFile];
            int bestSN = best.ScanNumber;
            int first = firstLastSN[rawFile][0];
            int last = firstLastSN[rawFile][1];
            ThermoRawFileScan current;
            double firstScanNum; //First scan in MS1 range
            double lastScanNum; //Last scan in MS1 range
            double firstTime = rawFile.ScanRange.FirstScanTime; //Earliest time in run
            double lastTime = rawFile.ScanRange.LastScanTime; //Latest time in run

            if (Form1.NEUCODE_12PLEX && Form1.MULTIINJECT)
            {
                fullScanList = new List<int>();
                foreach (PeptideSpectralMatch psm in allPSMs)
                {
                    current = rawFile[psm.ScanNumber].GetNextMsnScan(1).GetNextMsnScan(2);
                    fullScanList.Add(current.ScanNum);
                }
            }
            else
            {
                fullScanList = new List<int>();
                ThermoRawFileScan firstScan = rawFile[first].GetPreviousMsnScan(1);
                ThermoRawFileScan lastScan = rawFile[last].GetPreviousMsnScan(1);
                ThermoRawFileScan bestScan = rawFile[bestSN].GetPreviousMsnScan(1);

                //For peptides eluting across more than 1000 scans, create retention time window around best PSM
                if (last - first < 1000)
                {
                    rTRange = new Range(firstScan.ScanTime - rTWindow, lastScan.ScanTime + rTWindow);
                }
                else
                {
                    rTRange = new Range(bestScan.ScanTime - rTWindow, bestScan.ScanTime + rTWindow);
                    firstScan = bestScan;
                    lastScan = bestScan;
                }

                //Check to see if retention time window extends past the length of the run
                if (rTRange.MinValue < firstTime)
                {
                    rTRange.MinValue = firstTime; //Minimum time surpassed
                }

                if (rTRange.MaxValue > lastTime)
                {
                    rTRange.MaxValue = lastTime; //Maximum time surpassed
                }

                //Find the appropriate first and last full scan numbers
                current = firstScan;

                //For applications including 30K & 480K scans, skip over 30K MS1 scans
                if (Form1.NEUCODE && !Form1.CALCIUM)
                {
                    try
                    {
                        //Set first scan number
                        while (current.GetPreviousMsnScan(1).GetPreviousMsnScan(1) != null && current.GetPreviousMsnScan(1).GetPreviousMsnScan(1).ScanTime > rTRange.MinValue)
                        {
                            current = current.GetPreviousMsnScan(1).GetPreviousMsnScan(1);
                        }
                        firstScanNum = current.ScanNum;
                    }
                    catch (NullReferenceException)
                    {
                        Console.WriteLine("beginning of run exceeded");
                        firstScanNum = current.ScanNum;
                    }

                    //Set last scan number
                    current = lastScan;
                    try
                    {
                        while (current.GetNextMsnScan(1).GetNextMsnScan(1) != null && current.GetNextMsnScan(1).GetNextMsnScan(1).ScanTime < rTRange.MaxValue)
                        {
                            current = current.GetNextMsnScan(1).GetNextMsnScan(1);
                        }
                        lastScanNum = current.ScanNum;
                    }
                    catch (NullReferenceException)
                    {
                        Console.WriteLine("end of run exceeded");
                        lastScanNum = current.ScanNum;
                    }
                }
                //Otherwise, include every MS1 scan for quantification
                else
                {
                    //Set first scan number
                    while (current.GetPreviousMsnScan(1) != null && current.GetPreviousMsnScan(1).ScanTime > rTRange.MinValue)
                    {
                        current = current.GetPreviousMsnScan(1);
                    }
                    firstScanNum = current.ScanNum;

                    //Set last scan number
                    current = lastScan;
                    while (current.GetNextMsnScan(1) != null && current.GetNextMsnScan(1).ScanTime < rTRange.MaxValue)
                    {
                        current = current.GetNextMsnScan(1);
                    }
                    lastScanNum = current.ScanNum;
                }

                //Set the peptide's full scan range
                fullScanRange = new Range(firstScanNum, lastScanNum);

                //Extract scan numbers for list
                current = rawFile[(int)firstScanNum];

                while (current != null && current.ScanNum <= (int)fullScanRange.MaxValue)
                {
                    fullScanList.Add(current.ScanNum);
                    if (Form1.NEUCODE)
                    {
                        try
                        {
                            current = current.GetNextMsnScan(1).GetNextMsnScan(1);
                        }
                        catch (NullReferenceException)
                        {
                            Console.WriteLine("end of run exceeded");
                            break;
                        }
                    }
                    else
                    {
                        current = current.GetNextMsnScan(1);
                    }
                }
            }
        }

        /* Finds light and heavy peaks in the current scan, using a specified ppm window centered around the peptide's adjusted masses
         * Adds a pair to the peptide's list of all pairs for that raw file if enough (i.e., equal to or greater than the number of channels) peaks are found
         */
        public void findPeaks(ThermoRawFileScan current, ThermoRawFile rawFile, int charge)
        {
            peaks = new LabeledPeak[numChannels, numIsotopes];
            Pair pair;
            Tolerance tolerance;
            double SN = Form1.MINIMUMSN;
            int peaksNeeded;
            if (Form1.NOISEBANDCAP)
            {
                if (!Form1.NEUCODE)
                {
                    peaksNeeded = numChannels / 2;
                }
                else
                {
                    peaksNeeded = numIsotopologues / 2;
                }
            }
            else
            {
                if (!Form1.NEUCODE)
                {
                    peaksNeeded = numChannels;
                }
                else
                {
                    peaksNeeded = numIsotopologues;
                }
            }

            if (Form1.NEUCODE) //Use 5ppm tolerance for OMNE (480K resolution)
            {
                tolerance = NEUCODE;
            }
            else //Use 10ppm tolerance for SILAC (30K resolution)
            {
                tolerance = SILAC;
            }

            //Count the number of non-null peaks returned
            int peaksCount = 0;

            if (Form1.NEUCODE_ARG_PROLINECONVERSION)
            {
                double prolineConversion = 5.016775;
                LabeledPeak prolinePeak;
                for (int j = 0; j < numIsotopes; j++)
                {
                    peaks[0, j] = largestPeak(adjustedTheoMasses[0, j], current, tolerance, rawFile);
                    if (peaks[0, j] != null)
                    {
                        peaksCount++;
                        if (prolineCount > 0)
                        {
                            for (int p = 1; p <= prolineCount; p++)
                            {
                                prolinePeak = largestPeak(adjustedTheoMasses[0, j] + (p * prolineConversion), current, tolerance, rawFile);
                                if (prolinePeak != null)
                                {
                                    if (peaks[1, j] == null)
                                    {
                                        peaks[1, j] = prolinePeak;
                                    }
                                    else
                                    {
                                        peaks[1, j].Intensity += prolinePeak.Intensity;
                                    }
                                }
                            }
                        }
                        else
                        {
                            int fakeProlines = 2;
                            for (int p = 1; p <= fakeProlines; p++)
                            {
                                prolinePeak = largestPeak(adjustedTheoMasses[0, j] + (p * prolineConversion), current, tolerance, rawFile);
                                if (prolinePeak != null)
                                {
                                    if (peaks[1, j] == null)
                                    {
                                        peaks[1, j] = prolinePeak;
                                    }
                                    else
                                    {
                                        peaks[1, j].Intensity += prolinePeak.Intensity;
                                    }
                                }
                            }
                        }
                    }
                    else
                    {
                        peaks[0, j] = null;
                    }
                }

                double injectionTime = rawFile.GetInjectionTime(current.ScanNum);
                
                if (peaksCount > 0)
                {
                    for (int i = 0; i < 2; i++)
                    {
                        for (int j = 0; j < numIsotopes; j++)
                        {
                            if (peaks[i, j] != null)
                            {
                                peaks[i, j].dNL = peaks[i, j].Intensity * injectionTime;
                            }
                        }
                    }
                }

                //Create a new pair associated with the MS1 scan and peptide
                pair = new Pair(this, rawFile, current.ScanNum);
                pair.peaks = peaks;

                //Only add that pair to the quantification list if the number of non-null peaks returned is equal to or greater than the number of channels
                if (peaksCount >= 1)
                {
                    allHILACPairs[rawFile].Add(pair);
                }
            }
            /*else if (numChannels == 12)
            {
                List<LabeledPeak> peakList = new List<LabeledPeak>();
                for (int j = 0; j < numIsotopes; j++)
                {
                    for (int i = 0; i < numChannels; i++)
                    {
                        peaks[i, j] = largestPeak(adjustedTheoMasses[i, j], current, tolerance, rawFile);
                        if (peaks[i, j] != null)
                        {
                            peaksCount++;
                            peakList.Add(peaks[i, j]);
                        }
                    }
                }

                //Only add that pair to the quantification list if the number of non-null peaks returned is equal to or greater than the number of channels
                if (peaksCount >= peaksNeeded)
                {
                    //Create a new pair associated with the MS1 scan and peptide
                    pair = new Pair(this, rawFile, current.ScanNum);
                    pair.peaks = peaks; 
                    allHILACPairs[rawFile].Add(pair);
                }
                
                if (peaksCount == 0)
                {
                    Form1.ZEROPEAKS_COUNT++;
                }
                else if (peaksCount == 1)
                {
                    Form1.ONEPEAK_COUNT++;
                }

                
                else if (peaksCount == 2)
                {
                    Form1.TWOPEAKS_COUNT++;
                    peakList.Sort(LabeledPeak.sortMZAscending());
                    if (peakList[0].MZ != peakList[1].MZ)
                    {
                        spacing = ((peakList[1].MZ - peakList[0].MZ) * (double) charge) / (double) numLabels;
                        Form1.TWOPEAKS_SPACE.Add(spacing);
                    }
                    else
                    {  
                        Form1.TWOPEAKS_DUPLICATE++;
                    }
                }
                else if (peaksCount == 3)
                {
                    Form1.THREEPEAKS_COUNT++;
                    peakList.Sort(LabeledPeak.sortMZAscending());
                    if (peakList[0].MZ != peakList[1].MZ)
                    {
                        spacing = ((peakList[1].MZ - peakList[0].MZ) * (double)charge) / (double)numLabels;
                        Form1.THREEPEAKS_SPACE.Add(spacing);
                    }
                    else
                    {
                        Form1.THREEPEAKS_DUPLICATE++;
                    }
                    if (peakList[1].MZ != peakList[2].MZ)
                    {
                        spacing = ((peakList[2].MZ - peakList[1].MZ) * (double)charge) / (double)numLabels;
                        Form1.THREEPEAKS_SPACE.Add(spacing);
                    }
                    else
                    {
                        Form1.THREEPEAKS_DUPLICATE++;
                    }
                }
                else //peaksCount == 4
                {
                    Form1.FOURPEAKS_COUNT++;
                    peakList.Sort(LabeledPeak.sortMZAscending());
                    if (peakList[0].MZ != peakList[1].MZ)
                    {
                        spacing = ((peakList[1].MZ - peakList[0].MZ) * (double)charge) / (double)numLabels;
                        Form1.FOURPEAKS_SPACE.Add(spacing);
                    }
                    else
                    {
                        Form1.FOURPEAKS_DUPLICATE++;
                    }
                    if (peakList[1].MZ != peakList[2].MZ)
                    {
                        spacing = ((peakList[2].MZ - peakList[1].MZ) * (double)charge) / (double)numLabels;
                        Form1.FOURPEAKS_SPACE.Add(spacing);
                    }
                    else
                    {
                        Form1.FOURPEAKS_DUPLICATE++;
                    }
                    if (peakList[2].MZ != peakList[3].MZ)
                    {
                        spacing = ((peakList[3].MZ - peakList[2].MZ) * (double)charge) / (double)numLabels;
                        Form1.FOURPEAKS_SPACE.Add(spacing);
                    }
                    else
                    {
                        Form1.FOURPEAKS_DUPLICATE++;
                    }
                }
            }*/
            else if (Form1.NEUCODE_DUPLEX_LYS || Form1.NEUCODE_SIXPLEX_MTRAQ)
            {
                List<LabeledPeak> nonNullPeaks;
                List<int> duplicatePositions;
                List<int> missingChannel;
                int channelIndex;
                LabeledPeak peak;
                int overallPeaksCount = 0;
                double injectionTime = rawFile.GetInjectionTime(current.ScanNum);
                for (int c = 0; c < numClusters; c++) // Start cluster loop
                {
                    for (int j = 0; j < numIsotopes; j++) // Start isotope loop
                    {
                        nonNullPeaks = new List<LabeledPeak>();
                        duplicatePositions = new List<int>();
                        missingChannel = new List<int>();
                        channelIndex = c * numIsotopologues;

                        for (int i = 0; i < numIsotopologues; i++) // Start isotopologue loop
                        {
                            peak = largestPeak(adjustedTheoMasses[channelIndex + i, j], current, tolerance, rawFile);
                            if (peak != null)
                            {
                                nonNullPeaks.Add(peak);
                                peaks[channelIndex + i, j] = peak;
                            }
                            else
                            {
                                missingChannel.Add(i);
                            }
                        }

                        // Check for duplicate peaks
                        if (nonNullPeaks.Count > 1)
                        {
                            for (int m = 0; m < nonNullPeaks.Count - 1; m++)
                            {
                                if (nonNullPeaks[m].MZ == nonNullPeaks[m + 1].MZ)
                                {
                                    duplicatePositions.Add(m);
                                    duplicatePositions.Add(m + 1);
                                }
                            }
                        }

                        // Search for patterns if missing or duplicate peaks detected

                        if (nonNullPeaks.Count > 1 && (duplicatePositions.Count > 1 || missingChannel.Count > 0))
                        {
                            double lowerTolerance = tolerance.Value;
                            double upperTolerance = tolerance.Value;
                            double minMZ = Mass.MZFromMass(adjustedTheoMasses[channelIndex, j], charge);
                            double min = minMZ - Tolerance.GetThfromPPM(lowerTolerance, minMZ);
                            double maxMZ = Mass.MZFromMass(adjustedTheoMasses[channelIndex + 1, j], charge);
                            double max = maxMZ + Tolerance.GetThfromPPM(upperTolerance, maxMZ);
                            List<IPeak> peaksReTry = null;
                            List<LabeledPeak> sortedPeaks = null;
                            List<LabeledPeak> top2Peaks = null;
                            if (current.Spectrum.TryGetPeaks(out peaksReTry, min, max))
                            {
                                sortedPeaks = new List<LabeledPeak>();
                                foreach (IPeak peakReTry in peaksReTry)
                                {
                                    LabeledPeak newPeak = (LabeledPeak)peakReTry;
                                    if (newPeak.SN < SN)
                                    {
                                        newPeak = null;
                                    }
                                    if (newPeak != null && newPeak.Charge != 0 && newPeak.Charge != charge)
                                    {
                                        newPeak = null;
                                    }
                                    if (newPeak != null)
                                    {
                                        sortedPeaks.Add(newPeak);
                                    }
                                }
                                sortedPeaks.Sort(LabeledPeak.sortIntensityDescending());

                                // Peaks found for all 2 isotopologues
                                if (sortedPeaks.Count > 1)
                                {
                                    top2Peaks = new List<LabeledPeak>();
                                    top2Peaks.Add(sortedPeaks[0]);
                                    top2Peaks.Add(sortedPeaks[1]);
                                    top2Peaks.Sort(LabeledPeak.sortMZAscending());
                                    peaks[channelIndex, j] = top2Peaks[0];
                                    peaks[channelIndex + 1, j] = top2Peaks[1];
                                }
                                else if (sortedPeaks.Count == 1)
                                {
                                    LabeledPeak[] mappedPeaks = mapPeaks(sortedPeaks, channelIndex, channelIndex + 1, j, charge);
                                    for (int m = channelIndex; m < channelIndex + numIsotopologues; m++)
                                    {
                                        peaks[m, j] = mappedPeaks[m - channelIndex];
                                    }
                                }
                                else
                                {
                                    peaks[channelIndex, j] = null;
                                    peaks[channelIndex + 1, j] = null;
                                }
                            }
                        } // End pattern detection

                        peaksCount = 0;

                        // Count detected peaks
                        for (int i = channelIndex; i < channelIndex + numIsotopologues; i++)
                        {
                            if (peaks[i, j] != null)
                            {
                                peaks[i, j].dNL = peaks[i, j].Intensity * injectionTime;
                                peaksCount++;
                            }
                        }

                        if (peaksCount >= peaksNeeded)
                        {
                            overallPeaksCount++;
                        }
                        else
                        {
                            for (int i = channelIndex; i < channelIndex + numIsotopologues; i++)
                            {
                                peaks[i, j] = null;
                            }
                        }
                    } // End isotope loop
                } // End cluster loop

                if (overallPeaksCount > 0)
                {
                    pair = new Pair(this, rawFile, current.ScanNum);
                    pair.peaks = peaks;
                    allHILACPairs[rawFile].Add(pair);
                    quantified = true;
                }
            }
            else if (numIsotopologues == 4)
            {
                List<LabeledPeak> nonNullPeaks;
                List<int> duplicatePositions;
                List<int> missingChannel;
                int channelIndex;
                LabeledPeak peak;
                int overallPeaksCount = 0;
                double injectionTime = rawFile.GetInjectionTime(current.ScanNum);
                for (int c = 0; c < numClusters; c++) // Start cluster loop
                {
                    for (int j = 0; j < numIsotopes; j++) // Start isotope loop
                    {
                        nonNullPeaks = new List<LabeledPeak>();
                        duplicatePositions = new List<int>();
                        missingChannel = new List<int>();
                        channelIndex = c * numIsotopologues;

                        for (int i = 0; i < numIsotopologues; i++) // Start isotopologue loop
                        {
                            peak = largestPeak(adjustedTheoMasses[channelIndex + i, j], current, tolerance, rawFile);
                            if (peak != null)
                            {
                                nonNullPeaks.Add(peak);
                                peaks[channelIndex + i, j] = peak;
                            }
                            else
                            {
                                missingChannel.Add(i);
                            }
                        }

                        // Check for duplicate peaks
                        if (nonNullPeaks.Count > 1)
                        {
                            for (int m = 0; m < nonNullPeaks.Count - 1; m++)
                            {
                                if (nonNullPeaks[m].MZ == nonNullPeaks[m + 1].MZ)
                                {
                                    duplicatePositions.Add(m);
                                    duplicatePositions.Add(m + 1);
                                }
                            }
                        }

                        // Search for patterns if missing or duplicate peaks detected

                        if (duplicatePositions.Count > 1 || missingChannel.Count > 0)
                        {
                            double lowerTolerance = tolerance.Value;
                            double upperTolerance = tolerance.Value;
                            double minMZ = Mass.MZFromMass(adjustedTheoMasses[channelIndex, j], charge);
                            double min = minMZ - Tolerance.GetThfromPPM(lowerTolerance, minMZ);
                            double maxMZ = Mass.MZFromMass(adjustedTheoMasses[channelIndex + 3, j], charge);
                            double max = maxMZ + Tolerance.GetThfromPPM(upperTolerance, maxMZ);
                            List<IPeak> peaksReTry = null;
                            List<LabeledPeak> sortedPeaks = null;
                            List<LabeledPeak> top4Peaks = null;
                            if (current.Spectrum.TryGetPeaks(out peaksReTry, min, max))
                            {
                                sortedPeaks = new List<LabeledPeak>();
                                foreach (IPeak peakReTry in peaksReTry)
                                {
                                    LabeledPeak newPeak = (LabeledPeak)peakReTry;
                                    if (newPeak.SN < SN)
                                    {
                                        newPeak = null;
                                    }
                                    if (newPeak != null && newPeak.Charge != 0 && newPeak.Charge != charge)
                                    {
                                        newPeak = null;
                                    }
                                    if (newPeak != null)
                                    {
                                        sortedPeaks.Add(newPeak);
                                    }
                                }
                                sortedPeaks.Sort(LabeledPeak.sortIntensityDescending());

                                // Peaks found for all 4 isotopologues
                                if (sortedPeaks.Count > 3)
                                {
                                    top4Peaks = new List<LabeledPeak>();
                                    top4Peaks.Add(sortedPeaks[0]);
                                    top4Peaks.Add(sortedPeaks[1]);
                                    top4Peaks.Add(sortedPeaks[2]);
                                    top4Peaks.Add(sortedPeaks[3]);
                                    top4Peaks.Sort(LabeledPeak.sortMZAscending());
                                    peaks[channelIndex, j] = top4Peaks[0];
                                    peaks[channelIndex + 1, j] = top4Peaks[1];
                                    peaks[channelIndex + 2, j] = top4Peaks[2];
                                    peaks[channelIndex + 3, j] = top4Peaks[3];
                                }

                                // Only deal with fewer than 4 peaks when noise-band capping is enabled
                                else if (Form1.NOISEBANDCAP) 
                                {
                                    List<LabeledPeak> patternPeaksToMap = new List<LabeledPeak>();
                                    List<LabeledPeak> nonNullPeaksToMap = new List<LabeledPeak>();
                                    int duplicates = 0;
                                    foreach (LabeledPeak sortedPeak in sortedPeaks)
                                    {
                                        if (sortedPeak.SN >= 5.0)
                                        {
                                            patternPeaksToMap.Add(sortedPeak);
                                        }
                                    }
                                    foreach(LabeledPeak nonNullPeak in nonNullPeaks)
                                    {
                                        if (nonNullPeak.SN >= 5.0)
                                        {
                                            if (!nonNullPeaksToMap.Contains(nonNullPeak))
                                            {
                                                nonNullPeaksToMap.Add(nonNullPeak);
                                            }       
                                        }
                                    }
                                    if (patternPeaksToMap.Count > 1 && patternPeaksToMap.Count > nonNullPeaksToMap.Count)
                                    {
                                        patternPeaksToMap.Sort(LabeledPeak.sortMZAscending());
                                        LabeledPeak[] mappedPeaks = mapPeaks(patternPeaksToMap, channelIndex, channelIndex + (numIsotopologues - 1), j, charge);
                                        for (int m = channelIndex; m < channelIndex + numIsotopologues; m++)
                                        {
                                            peaks[m, j] = mappedPeaks[m - channelIndex];
                                        }
                                    }
                                    else if (nonNullPeaksToMap.Count > 1)
                                    {
                                        nonNullPeaksToMap.Sort(LabeledPeak.sortMZAscending());
                                        LabeledPeak[] mappedPeaks = mapPeaks(nonNullPeaksToMap, channelIndex, channelIndex + (numIsotopologues - 1), j, charge);
                                        for (int m = channelIndex; m < channelIndex + numIsotopologues; m++)
                                        {
                                            peaks[m, j] = mappedPeaks[m - channelIndex];
                                        }
                                    }
                                    else
                                    {
                                        peaks[channelIndex, j] = null;
                                        peaks[channelIndex + 1, j] = null;
                                        peaks[channelIndex + 2, j] = null;
                                        peaks[channelIndex + 3, j] = null;
                                    }
                                }
                                else
                                {
                                    peaks[channelIndex, j] = null;
                                    peaks[channelIndex + 1, j] = null;
                                    peaks[channelIndex + 2, j] = null;
                                    peaks[channelIndex + 3, j] = null;
                                }
                            }
                        } // End pattern detection

                        peaksCount = 0;

                        // Count detected peaks
                        for (int i = channelIndex; i < channelIndex + numIsotopologues; i++)
                        {
                            if (peaks[i, j] != null)
                            {
                                peaks[i, j].dNL = peaks[i, j].Intensity * injectionTime;
                                peaksCount++;
                            }
                        }

                        if (peaksCount >= peaksNeeded)
                        {
                            overallPeaksCount++;
                        }
                        else
                        {
                            for (int i = channelIndex; i < channelIndex + numIsotopologues; i++)
                            {
                                peaks[i, j] = null;
                            }
                        }
                    } // End isotope loop
                } // End cluster loop

                if (overallPeaksCount > 0)
                {
                    pair = new Pair(this, rawFile, current.ScanNum);
                    pair.peaks = peaks;
                    allHILACPairs[rawFile].Add(pair);
                    quantified = true;
                }
            }

            else if (Form1.NEUCODE_DUPLEX_LYS1 || Form1.NEUCODE_DUPLEX_CARBAMYL)
            {
                for (int j = 0; j < numIsotopes; j++)
                {
                    for (int i = 0; i < numChannels; i += 2)
                    {
                        List<LabeledPeak> nonNullPeaks = new List<LabeledPeak>();
                        List<int> duplicatePositions = new List<int>();
                        peaks[i, j] = largestPeak(adjustedTheoMasses[i, j], current, tolerance, rawFile);
                        peaks[i + 1, j] = largestPeak(adjustedTheoMasses[i + 1, j], current, tolerance, rawFile);

                        if (peaks[i, j] != null)
                        {
                            nonNullPeaks.Add(peaks[i, j]);
                            peaksCount++;
                        }
                        if (peaks[i + 1, j] != null)
                        {
                            nonNullPeaks.Add(peaks[i + 1, j]);
                            peaksCount++;
                        }

                        if (nonNullPeaks.Count > 1)
                        {
                            for (int m = 0; m < nonNullPeaks.Count - 1; m++)
                            {
                                if (nonNullPeaks[m].MZ == nonNullPeaks[m + 1].MZ)
                                {
                                    duplicatePositions.Add(m);
                                    duplicatePositions.Add(m + 1);
                                }
                            }
                        }

                        if (duplicatePositions.Count > 1)
                        {
                            double minMZ = Mass.MZFromMass(adjustedTheoMasses[i, j], charge);
                            double min = minMZ - Tolerance.GetThfromPPM(10.0, minMZ);
                            double maxMZ = Mass.MZFromMass(adjustedTheoMasses[i + 1, j], charge);
                            double max = maxMZ + Tolerance.GetThfromPPM(10.0, maxMZ);
                            List<IPeak> peaksReTry = null;
                            List<LabeledPeak> sortedPeaks = null;
                            List<LabeledPeak> top2Peaks = null;
                            if (current.Spectrum.TryGetPeaks(out peaksReTry, min, max))
                            {
                                sortedPeaks = new List<LabeledPeak>();
                                foreach (IPeak peakReTry in peaksReTry)
                                {
                                    LabeledPeak newPeak = (LabeledPeak)peakReTry;
                                    if (newPeak.SN < Form1.MINIMUMSN)
                                    {
                                        newPeak = null;
                                    }
                                    if (newPeak != null && newPeak.Charge != 0 && newPeak.Charge != charge)
                                    {
                                        newPeak = null;
                                    }
                                    if (newPeak != null)
                                    {
                                        sortedPeaks.Add(newPeak);
                                    }
                                }
                                sortedPeaks.Sort(LabeledPeak.sortIntensityDescending());

                                if (sortedPeaks.Count > 1)
                                {
                                    top2Peaks = new List<LabeledPeak>();
                                    top2Peaks.Add(sortedPeaks[0]);
                                    top2Peaks.Add(sortedPeaks[1]);
                                    top2Peaks.Sort(LabeledPeak.sortMZAscending());
                                    peaks[i, j] = top2Peaks[0];
                                    peaks[i + 1, j] = top2Peaks[1];
                                    peaksCount = peaksCount - nonNullPeaks.Count + 2;
                                }
                                else
                                {
                                    LabeledPeak[] mappedPeaks = mapPeaks(nonNullPeaks, 0, 1, j, charge);
                                    for (int m = 0; m < mappedPeaks.Length; m++)
                                    {
                                        peaks[m, j] = mappedPeaks[m];
                                    }
                                    peaksCount = peaksCount - (duplicatePositions.Count / 2);
                                }
                            }
                        }
                    }

                    //Create a new pair associated with the MS1 scan and peptide
                    pair = new Pair(this, rawFile, current.ScanNum);
                    pair.peaks = peaks;

                    //Only add that pair to the quantification list if the number of non-null peaks returned is equal to or greater than the number of channels
                    if (peaksCount >= peaksNeeded)
                    {
                        double injectionTime = rawFile.GetInjectionTime(current.ScanNum);
                        for (int i = 0; i < numChannels; i++)
                        {
                            for (int k = 0; k < numIsotopes; k++)
                            {
                                if (pair.peaks[i, k] != null)
                                {
                                    pair.peaks[i, k].dNL = pair.peaks[i, k].Intensity * injectionTime;
                                }
                            }
                        }

                        allHILACPairs[rawFile].Add(pair);
                    }
                }
            }
            else
            {
                //double signalNoise = 0;
                for (int j = 0; j < numIsotopes; j++)
                {
                    for (int i = 0; i < numChannels; i++)
                    {
                        peaks[i, j] = largestPeak(adjustedTheoMasses[i, j], current, tolerance, rawFile);
                        if (peaks[i, j] != null)
                        {
                            peaksCount++;
                            //signalNoise += peaks[i,j].SN;
                        }
                        else
                        {
                            peaks[i, j] = null;
                        }
                    }
                }
                /*if (peaksCount > 0)
                {
                    List<double> SNList;
                    if (Form1.PEAKSNPAIRS.TryGetValue(peaksCount, out SNList))
                    {
                        SNList.Add(signalNoise / (double)peaksCount);
                    }
                    else
                    {
                        SNList = new List<double>();
                        SNList.Add(signalNoise / (double)peaksCount);
                        Form1.PEAKSNPAIRS.Add(peaksCount, SNList);
                    }
                }*/


                //Create a new pair associated with the MS1 scan and peptide
                pair = new Pair(this, rawFile, current.ScanNum);
                pair.peaks = peaks;

                //Only add that pair to the quantification list if the number of non-null peaks returned is equal to or greater than the number of channels
                if (peaksCount >= numChannels)
                {
                    allHILACPairs[rawFile].Add(pair);
                }
            }
        }

        public LabeledPeak[] mapPeaks(List<LabeledPeak> allPeaks, int channelStart, int channelEnd, int isotope, int charge)
        {
            int numChannels = channelEnd - channelStart + 1;
            LabeledPeak[] mappedPeaks = new LabeledPeak[numChannels];
            List<LabeledPeak> peaks = new List<LabeledPeak>();

            //First remove any duplicate peaks
            foreach (LabeledPeak peak in allPeaks)
            {
                if (peaks.Count == 0 || !peaks.Contains(peak))
                {
                    peaks.Add(peak);
                }
            }

            peaks.Sort(LabeledPeak.sortMZAscending());
            
            int possibleCombinations;
            if (peaks.Count == 2)
            {
                possibleCombinations = 6;
            }
            else if (peaks.Count == 3)
            {
                possibleCombinations = 4;
            }
            else
            {
                possibleCombinations = 2;
            }
            KeyValuePair<double, LabeledPeak>[,] errorCheck = new KeyValuePair<double, LabeledPeak>[possibleCombinations, numChannels];
            
            //Then calculate all the errors
            LabeledPeak currentPeak;
            double neutralMass;
            double error1;
            double error2;
            double error3;
            double error4;
            double overallError;
            double minError;
            int minErrorIndex;
            for (int i = 0; i < peaks.Count; i++) //Peaks loop
            {
                currentPeak = peaks[i];

                if (Form1.NEUCODE_DUPLEX_LYS1 || Form1.NEUCODE_DUPLEX_CARBAMYL || Form1.NEUCODE_SIXPLEX_MTRAQ || Form1.NEUCODE_DUPLEX_LYS)
                {
                    neutralMass = Mass.MassFromMz(currentPeak.MZ, charge);

                    if (Form1.NUMCHANNELS == 2)
                    {
                        error1 = Tolerance.GetError(neutralMass, adjustedTheoMasses[0, isotope], ToleranceType.PPM);
                        error2 = Tolerance.GetError(neutralMass, adjustedTheoMasses[1, isotope], ToleranceType.PPM);
                        if (Math.Abs(error1) < Math.Abs(error2))
                        {
                            mappedPeaks[0] = currentPeak;
                        }
                        else
                        {
                            mappedPeaks[1] = currentPeak;
                        }
                    }
                    else
                    {
                        error1 = Tolerance.GetError(neutralMass, adjustedTheoMasses[channelStart, isotope], ToleranceType.PPM);
                        error2 = Tolerance.GetError(neutralMass, adjustedTheoMasses[channelEnd, isotope], ToleranceType.PPM);
                        if (Math.Abs(error1) < Math.Abs(error2))
                        {
                            mappedPeaks[0] = currentPeak;
                        }
                        else
                        {
                            mappedPeaks[1] = currentPeak;
                        }
                    }
                    return mappedPeaks;
                }
                else
                {
                    neutralMass = Mass.MassFromMz(currentPeak.MZ, charge);
                    error1 = Tolerance.GetError(neutralMass, adjustedTheoMasses[channelStart, isotope], ToleranceType.PPM);
                    error2 = Tolerance.GetError(neutralMass, adjustedTheoMasses[channelStart + 1, isotope], ToleranceType.PPM);
                    error3 = Tolerance.GetError(neutralMass, adjustedTheoMasses[channelStart + 2, isotope], ToleranceType.PPM);
                    error4 = Tolerance.GetError(neutralMass, adjustedTheoMasses[channelStart + 3, isotope], ToleranceType.PPM);

                    if (peaks.Count == 2)
                    {
                        if (i == 0)
                        {
                            errorCheck[0, 0] = new KeyValuePair<double, LabeledPeak>(Math.Abs(error1), currentPeak);
                            errorCheck[1, 0] = new KeyValuePair<double, LabeledPeak>(Math.Abs(error1), currentPeak);
                            errorCheck[2, 0] = new KeyValuePair<double, LabeledPeak>(Math.Abs(error1), currentPeak);
                            errorCheck[3, 1] = new KeyValuePair<double, LabeledPeak>(Math.Abs(error2), currentPeak);
                            errorCheck[4, 1] = new KeyValuePair<double, LabeledPeak>(Math.Abs(error2), currentPeak);
                            errorCheck[5, 2] = new KeyValuePair<double, LabeledPeak>(Math.Abs(error3), currentPeak);
                        }
                        else
                        {
                            errorCheck[0, 1] = new KeyValuePair<double, LabeledPeak>(Math.Abs(error2), currentPeak);
                            errorCheck[1, 2] = new KeyValuePair<double, LabeledPeak>(Math.Abs(error3), currentPeak);
                            errorCheck[2, 3] = new KeyValuePair<double, LabeledPeak>(Math.Abs(error4), currentPeak);
                            errorCheck[3, 2] = new KeyValuePair<double, LabeledPeak>(Math.Abs(error3), currentPeak);
                            errorCheck[4, 3] = new KeyValuePair<double, LabeledPeak>(Math.Abs(error4), currentPeak);
                            errorCheck[5, 3] = new KeyValuePair<double, LabeledPeak>(Math.Abs(error4), currentPeak);
                        }
                    }

                    if (peaks.Count == 3)
                    {
                        if (i == 0)
                        {
                            errorCheck[0, 0] = new KeyValuePair<double, LabeledPeak>(Math.Abs(error1), currentPeak);
                            errorCheck[1, 0] = new KeyValuePair<double, LabeledPeak>(Math.Abs(error1), currentPeak);
                            errorCheck[2, 0] = new KeyValuePair<double, LabeledPeak>(Math.Abs(error1), currentPeak);
                            errorCheck[3, 1] = new KeyValuePair<double, LabeledPeak>(Math.Abs(error2), currentPeak);
                        }
                        else if (i == 1)
                        {
                            errorCheck[0, 1] = new KeyValuePair<double, LabeledPeak>(Math.Abs(error2), currentPeak);
                            errorCheck[1, 1] = new KeyValuePair<double, LabeledPeak>(Math.Abs(error2), currentPeak);
                            errorCheck[2, 2] = new KeyValuePair<double, LabeledPeak>(Math.Abs(error3), currentPeak);
                            errorCheck[3, 2] = new KeyValuePair<double, LabeledPeak>(Math.Abs(error3), currentPeak);
                        }
                        else
                        {
                            errorCheck[0, 2] = new KeyValuePair<double, LabeledPeak>(Math.Abs(error3), currentPeak);
                            errorCheck[1, 3] = new KeyValuePair<double, LabeledPeak>(Math.Abs(error4), currentPeak);
                            errorCheck[2, 3] = new KeyValuePair<double, LabeledPeak>(Math.Abs(error4), currentPeak);
                            errorCheck[3, 3] = new KeyValuePair<double, LabeledPeak>(Math.Abs(error4), currentPeak);
                        }
                    }
                }                
            }

            // Find combination that leads to lowest overall error
            minError = 100;
            minErrorIndex = -1;
            for (int j = 0; j < possibleCombinations; j++)
            {
                overallError = 0;
                for (int k = 0; k < numChannels; k++)
                {
                    overallError += errorCheck[j, k].Key;
                }

                if (overallError < minError)
                {
                    minError = overallError;
                    minErrorIndex = j;
                }
            }

            // Set mappedPeaks to combination with lowest overall error
            for (int n = 0; n < numChannels; n++)
            {
                mappedPeaks[n] = errorCheck[minErrorIndex, n].Value;
            }

            return mappedPeaks;
        }

        public void checkPeaks(ThermoRawFileScan current, ThermoRawFile rawFile, LabeledPeak[,] peaks, int scanNumber, Pair pair)
        {
            int numIsotopes = Form1.NUMISOTOPES;
            double injectionTime = rawFile.GetInjectionTime(current.ScanNum);
            double noise;
            List<double> noisePeaks = new List<double>();

            //Calculate the average noise
            for (int i = 0; i < numChannels; i++)
            {
                for (int j = 0; j < numIsotopes; j++)
                {
                    if (peaks[i,j] != null)
                    {
                        noisePeaks.Add(peaks[i,j].Noise * injectionTime);
                    }
                }
            }
            noise = noisePeaks.Sum() / (double) noisePeaks.Count;

            //Apply noise to missing channels
            for (int i = 0; i < numChannels; i++)
            {
                for (int j = 0; j < numIsotopes; j++)
                {
                    if (peaks[i, j] == null)
                    {
                        peaks[i, j] = new LabeledPeak();
                        peaks[i, j].Intensity = noise;
                        peaks[i, j].MZ = 0;
                    }
                }
            }
        }

        public void checkPairCoalescence(ThermoRawFile rawFile)
        {
            PeptideSpectralMatch best = PSMs[rawFile][0];
            int scanNumber = best.ScanNumber;
            ThermoRawFileScan current;

            List<Pair> pairs;
            Pair pair;
            LabeledPeak peak1;
            LabeledPeak peak2;
            LabeledPeak peak3;
            LabeledPeak peak4;
            IsotopePair isotopePair;
            List<IsotopePair> lowerIntensity;
            List<IsotopePair> higherIntensity;
            double lowerFrequency;
            double higherFrequency;

            try
            {
                if (Form1.MULTIINJECT)
                {
                    current = rawFile[scanNumber].GetNextMsnScan(1).GetNextMsnScan(2);
                }
                else
                {
                    current = rawFile[scanNumber].GetPrecursorRawScan();
                }

                allHILACPairs.TryGetValue(rawFile, out pairs);
                for (int p = 0; p < pairs.Count(); p++)
                {
                    pair = pairs[p];
                    for (int c = 0; c < numChannels; c += numIsotopologues)
                    {
                        for (int j = 0; j < numIsotopes; j++)
                        {
                            if (numIsotopologues == 2)
                            {
                                peak1 = pair.peaks[c, j];
                                peak2 = pair.peaks[c + 1, j];

                                if (peak1 != null && peak2 != null)
                                {
                                    // Do nothing for complete pairs
                                }
                                else if (peak1 == null && peak2 == null)
                                {
                                    // Do nothing for null pairs
                                }
                                else if (maximumIntensity <= System.Math.Log10(Form1.MAXIMUMDNL))
                                {
                                    // Do nothing for peptides below the threshold
                                }
                                else
                                {
                                    LabeledPeak singlePeak;
                                    LabeledPeak newPeak1;
                                    LabeledPeak newPeak2;
                                    double singlePeakIntensity;

                                    //Find the single peak
                                    if (peak1 == null)
                                    {
                                        singlePeak = peak2;
                                        singlePeakIntensity = peak2.dNL;
                                    }
                                    else
                                    {
                                        singlePeak = peak1;
                                        singlePeakIntensity = peak1.dNL;
                                    }

                                    //Sort the peptide's other pairs (either both or one channel present) based on their intensities relative to the single peak
                                    lowerIntensity = new List<IsotopePair>();
                                    higherIntensity = new List<IsotopePair>();
                                    Pair otherPair;
                                    for (int a = 0; a < pairs.Count(); a++)
                                    {
                                        otherPair = pairs[a];
                                        for (int h = 0; h < numChannels; h += numIsotopologues)
                                        {
                                            for (int s = 0; s < numIsotopes; s++)
                                            {
                                                newPeak1 = otherPair.peaks[h, s];
                                                newPeak2 = otherPair.peaks[h + 1, s];

                                                //Both channels present
                                                if (newPeak1 != null && newPeak2 != null)
                                                {
                                                    isotopePair = new IsotopePair(otherPair, s, newPeak1.MZ, newPeak1.dNL, newPeak2.MZ, newPeak2.dNL);
                                                    if (isotopePair.intensity <= singlePeak.dNL)
                                                    {
                                                        lowerIntensity.Add(isotopePair);
                                                    }
                                                    else
                                                    {
                                                        higherIntensity.Add(isotopePair);
                                                    }
                                                }
                                                //One channel present
                                                else
                                                {
                                                    if (newPeak1 != null)
                                                    {
                                                        isotopePair = new IsotopePair(otherPair, s, newPeak1.MZ, newPeak1.dNL, 0, 0);
                                                        if (isotopePair.intensity <= singlePeak.dNL)
                                                        {
                                                            lowerIntensity.Add(isotopePair);
                                                        }
                                                        else
                                                        {
                                                            higherIntensity.Add(isotopePair);
                                                        }
                                                    }
                                                    if (newPeak2 != null)
                                                    {
                                                        isotopePair = new IsotopePair(otherPair, s, newPeak2.MZ, newPeak2.dNL, 0, 0);
                                                        if (isotopePair.intensity <= singlePeak.dNL)
                                                        {
                                                            lowerIntensity.Add(isotopePair);
                                                        }
                                                        else
                                                        {
                                                            higherIntensity.Add(isotopePair);
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }

                                    //Calculate the missing channel frequencies for pairs both lower and higher in intensity than the single peak
                                    lowerFrequency = 0;
                                    int lowerCount = 0;
                                    for (int l = 0; l < lowerIntensity.Count(); l++)
                                    {
                                        if (lowerIntensity[l].missingChannel)
                                        {
                                            lowerCount++;
                                        }
                                    }
                                    lowerFrequency = (double)lowerCount / (double)lowerIntensity.Count();

                                    higherFrequency = 0;
                                    int higherCount = 0;
                                    for (int l = 0; l < higherIntensity.Count(); l++)
                                    {
                                        if (higherIntensity[l].missingChannel)
                                        {
                                            higherCount++;
                                        }
                                    }
                                    higherFrequency = (double)higherCount / (double)higherIntensity.Count();

                                    //A peak is deemed coalesced if the frequency of a missing channel is at least 1.5-fold greater for the more intense pairs than the less intense (each category must have at least 2 contributing pairs)
                                    if (lowerIntensity.Count() > 1 && higherIntensity.Count() > 1 && higherFrequency / lowerFrequency > 1.5)
                                    {
                                        if (coalescedPeakIntensities == null)
                                        {
                                            coalescedPeakIntensities = new List<double>();
                                        }
                                        else
                                        {
                                            coalescedPeakIntensities.Add(singlePeakIntensity);
                                        }

                                        //Set single peak to null
                                        if (peak1 == null)
                                        {
                                            pair.peaks[c + 1, j] = null;
                                        }
                                        if (peak2 == null)
                                        {
                                            pair.peaks[c, j] = null;
                                        }
                                    }
                                }
                            }

                            if (numIsotopologues == 4)
                            {
                                peak1 = pair.peaks[c, j];
                                peak2 = pair.peaks[c + 1, j];
                                peak3 = pair.peaks[c + 2, j];
                                peak4 = pair.peaks[c + 3, j];

                                if (peak1 != null && peak2 != null && peak3 != null && peak4 != null)
                                {
                                    // Do nothing for complete pairs
                                }
                                else if (peak1 == null && peak2 == null && peak3 == null && peak4 == null)
                                {
                                    // Do nothing for null pairs
                                }
                                else if (maximumIntensity <= System.Math.Log10(Form1.MAXIMUMDNL))
                                {
                                    // Do nothing for peptides below the threshold
                                }
                                else // Should have 2 or three peaks
                                {
                                    LabeledPeak newPeak1;
                                    LabeledPeak newPeak2;
                                    LabeledPeak newPeak3;
                                    LabeledPeak newPeak4;
                                    double presentPeaksIntensity = 0;

                                    // Calculate the maximum intensity of all present peaks
                                    if (peak1 != null)
                                    {
                                        if (peak1.dNL > presentPeaksIntensity)
                                        {
                                            presentPeaksIntensity = peak1.dNL;
                                        }
                                    }
                                    if (peak2 != null)
                                    {
                                        if (peak2.dNL > presentPeaksIntensity)
                                        {
                                            presentPeaksIntensity = peak2.dNL;
                                        }
                                    }
                                    if (peak3 != null)
                                    {
                                        if (peak3.dNL > presentPeaksIntensity)
                                        {
                                            presentPeaksIntensity = peak3.dNL;
                                        }
                                    }
                                    if (peak4 != null)
                                    {
                                        if (peak4.dNL > presentPeaksIntensity)
                                        {
                                            presentPeaksIntensity = peak4.dNL;
                                        }
                                    }

                                    //Sort the peptide's other pairs (either both or one channel present) based on their intensities relative to the single peak
                                    lowerIntensity = new List<IsotopePair>();
                                    higherIntensity = new List<IsotopePair>();
                                    Pair otherPair;
                                    for (int a = 0; a < pairs.Count(); a++)
                                    {
                                        otherPair = pairs[a];
                                        for (int h = 0; h < numChannels; h += numIsotopologues)
                                        {
                                            for (int s = 0; s < numIsotopes; s++)
                                            {
                                                newPeak1 = otherPair.peaks[h, s];
                                                newPeak2 = otherPair.peaks[h + 1, s];
                                                newPeak3 = otherPair.peaks[h + 2, s];
                                                newPeak4 = otherPair.peaks[h + 3, s];

                                                //Both channels present
                                                if (newPeak1 != null && newPeak2 != null && newPeak3 != null && newPeak4 != null)
                                                {
                                                    isotopePair = new IsotopePair(otherPair, s, newPeak1.MZ, newPeak1.dNL, newPeak2.MZ, newPeak2.dNL, newPeak3.MZ, newPeak3.dNL, newPeak4.MZ, newPeak4.dNL);
                                                    if (isotopePair.intensity <= presentPeaksIntensity)
                                                    {
                                                        lowerIntensity.Add(isotopePair);
                                                    }
                                                    else
                                                    {
                                                        higherIntensity.Add(isotopePair);
                                                    }
                                                }
                                                else if (newPeak1 == null && newPeak2 == null && newPeak3 == null && newPeak4 == null)
                                                {
                                                    // Do nothing for incomplete pairs
                                                }
                                                else
                                                {
                                                    LabeledPeak blank = new LabeledPeak();
                                                    blank.MZ = 0;
                                                    blank.dNL = 0;

                                                    if (newPeak1 == null)
                                                    {
                                                        newPeak1 = blank;
                                                    }
                                                    if (newPeak2 == null)
                                                    {
                                                        newPeak2 = blank;
                                                    }
                                                    if (newPeak3 == null)
                                                    {
                                                        newPeak3 = blank;
                                                    }
                                                    if (newPeak4 == null)
                                                    {
                                                        newPeak4 = blank;
                                                    }

                                                    isotopePair = new IsotopePair(otherPair, s, newPeak1.MZ, newPeak1.dNL, newPeak2.MZ, newPeak2.dNL, newPeak3.MZ, newPeak3.dNL, newPeak4.MZ, newPeak4.dNL);
                                                    if (isotopePair.intensity <= presentPeaksIntensity)
                                                    {
                                                        lowerIntensity.Add(isotopePair);
                                                    }
                                                    else
                                                    {
                                                        higherIntensity.Add(isotopePair);
                                                    }
                                                }
                                            }
                                        }
                                    }

                                    //Calculate the missing channel frequencies for pairs both lower and higher in intensity than the single peak
                                    lowerFrequency = 0;
                                    int lowerCount = 0;
                                    for (int l = 0; l < lowerIntensity.Count(); l++)
                                    {
                                        if (lowerIntensity[l].missingChannel)
                                        {
                                            lowerCount++;
                                        }
                                    }
                                    lowerFrequency = (double)lowerCount / (double)lowerIntensity.Count();

                                    higherFrequency = 0;
                                    int higherCount = 0;
                                    for (int l = 0; l < higherIntensity.Count(); l++)
                                    {
                                        if (higherIntensity[l].missingChannel)
                                        {
                                            higherCount++;
                                        }
                                    }
                                    higherFrequency = (double)higherCount / (double)higherIntensity.Count();

                                    //A peak is deemed coalesced if the frequency of a missing channel is at least 1.5-fold greater for the more intense pairs than the less intense (each category must have at least 2 contributing pairs)
                                    if (lowerIntensity.Count() > 1 && higherIntensity.Count() > 1 && higherFrequency / lowerFrequency > 1.5)
                                    {
                                        if (coalescedPeakIntensities == null)
                                        {
                                            coalescedPeakIntensities = new List<double>();
                                        }
                                        else
                                        {
                                            coalescedPeakIntensities.Add(presentPeaksIntensity);
                                        }

                                        //Set single peaks to null
                                        pair.peaks[c, j] = null;
                                        pair.peaks[c + 1, j] = null;
                                        pair.peaks[c + 2, j] = null;
                                        pair.peaks[c + 3, j] = null;
                                    }
                                }
                            }
                        } // End channel loop
                    } // End isotope loop
                } // End pair loop
            }
            catch (Exception)
            { }
        }

        public void checkPairSpacing(ThermoRawFile rawFile, List<Spacing> spacings)
        {
            PeptideSpectralMatch best = PSMs[rawFile][0];
            int charge = best.Charge;
            List<Pair> pairs;
            Pair pair;
            LabeledPeak peak1;
            LabeledPeak peak2;
            LabeledPeak peak3;
            LabeledPeak peak4;
            double spacing1;
            double spacing2;
            double spacing3;
            double peak1NeutralMass;
            double peak2NeutralMass;
            double peak3NeutralMass;
            double peak4NeutralMass;
            double lysineVersion;

            if (Form1.NEUCODE)
            {
                if (numChannels == 4 || numChannels == 12)
                {
                    lysineVersion = NEUCODE12;
                }
                else
                {
                    lysineVersion = NEUCODELYS;
                }
            }
            else
            {
                lysineVersion = heavyK.MonoisotopicMass;
            }

            /*if (Form1.NEUCODE_SIXPLEX_ARG && arginineCount == 0)
            {
                numChannels = 2;
            }*/

            allHILACPairs.TryGetValue(rawFile, out pairs);
            if (numIsotopologues == 4)
            {
                for (int p = 0; p < pairs.Count(); p++)
                {
                    pair = pairs[p];
                    int peakCount;
                    int overallPeakCount = 0;

                    for (int c = 0; c < numClusters; c++)
                    {
                        for (int j = 0; j < numIsotopes; j ++)
                        {
                            peakCount = 0;
                            int channelIndex = c * numIsotopologues;
                            for (int i = c; i < c + numIsotopologues; i += numIsotopologues)
                            { 
                                peak1 = pair.peaks[channelIndex, j];
                                peak2 = pair.peaks[channelIndex + 1, j];
                                peak3 = pair.peaks[channelIndex + 2, j];
                                peak4 = pair.peaks[channelIndex + 3, j];
                                if (peak1 != null)
                                {
                                    peakCount++;
                                }
                                if (peak2 != null)
                                {
                                    peakCount++;
                                }
                                if (peak3 != null)
                                {
                                    peakCount++;
                                }
                                if (peak4 != null)
                                {
                                    peakCount++;
                                }

                                //First consider complete pairs
                                if (peakCount == 4)
                                {
                                    peak1NeutralMass = Mass.MassFromMz(peak1.MZ, charge);
                                    peak2NeutralMass = Mass.MassFromMz(peak2.MZ, charge);
                                    peak3NeutralMass = Mass.MassFromMz(peak3.MZ, charge);
                                    peak4NeutralMass = Mass.MassFromMz(peak4.MZ, charge);

                                    spacing1 = peak4NeutralMass - peak3NeutralMass;
                                    spacing2 = peak3NeutralMass - peak2NeutralMass;
                                    spacing3 = peak2NeutralMass - peak1NeutralMass;

                                    if (Form1.NOISEBANDCAP)
                                    {
                                        if (spacing1 < spacingRange.MinValue || spacing1 > spacingRange.MaxValue)
                                        {
                                            pair.peaks[channelIndex + 3, j] = null;
                                        }
                                        if (spacing2 < spacingRange.MinValue || spacing2 > spacingRange.MaxValue)
                                        {
                                            pair.peaks[channelIndex + 2, j] = null;
                                        }
                                        if (spacing3 < spacingRange.MinValue || spacing3 > spacingRange.MaxValue)
                                        {
                                            pair.peaks[channelIndex, j] = null;
                                            pair.peaks[channelIndex + 1, j] = null;
                                        }
                                    }
                                    else
                                    {
                                        if (spacing1 < spacingRange.MinValue || spacing1 > spacingRange.MaxValue || spacing1 < spacingRange.MinValue || spacing1 > spacingRange.MaxValue || spacing3 < spacingRange.MinValue || spacing3 > spacingRange.MaxValue)
                                        {
                                            pair.peaks[channelIndex, j] = null;
                                            pair.peaks[channelIndex + 1, j] = null;
                                            pair.peaks[channelIndex + 2, j] = null;
                                            pair.peaks[channelIndex + 3, j] = null;
                                        }
                                    }
                                }
                                else if (peakCount == 3)
                                {
                                    if (Form1.NOISEBANDCAP)
                                    {
                                        if (peak1 != null)
                                        {
                                            peak1NeutralMass = Mass.MassFromMz(peak1.MZ, charge);
                                            if (peak2 != null && peak3 != null)
                                            {
                                                peak2NeutralMass = Mass.MassFromMz(peak2.MZ, charge);
                                                peak3NeutralMass = Mass.MassFromMz(peak3.MZ, charge);
                                                spacing1 = peak3NeutralMass - peak2NeutralMass;
                                                spacing2 = peak2NeutralMass - peak1NeutralMass;
                                                if (spacing1 > spacingRange.MaxValue || spacing1 < spacingRange.MinValue)
                                                {
                                                    pair.peaks[channelIndex + 2, j] = null;
                                                }
                                                if (spacing2 > spacingRange.MaxValue || spacing2 < spacingRange.MinValue)
                                                {
                                                    pair.peaks[channelIndex + 1, j] = null;
                                                    pair.peaks[channelIndex, j] = null;
                                                }
                                            }
                                            else if (peak3 != null && peak4 != null)
                                            {
                                                peak3NeutralMass = Mass.MassFromMz(peak3.MZ, charge);
                                                peak4NeutralMass = Mass.MassFromMz(peak4.MZ, charge);
                                                spacing1 = peak4NeutralMass - peak3NeutralMass;
                                                spacing2 = (peak3NeutralMass - peak1NeutralMass) / 2.0;
                                                if (spacing1 > spacingRange.MaxValue || spacing1 < spacingRange.MinValue)
                                                {
                                                    pair.peaks[channelIndex + 3, j] = null;
                                                }
                                                if (spacing2 > spacingRange.MaxValue || spacing2 < spacingRange.MinValue)
                                                {
                                                    pair.peaks[channelIndex + 2, j] = null;
                                                    pair.peaks[channelIndex, j] = null;
                                                }
                                            }
                                            else
                                            {
                                                peak2NeutralMass = Mass.MassFromMz(peak2.MZ, charge);
                                                peak4NeutralMass = Mass.MassFromMz(peak4.MZ, charge);
                                                spacing1 = (peak4NeutralMass - peak2NeutralMass) / 2.0;
                                                spacing2 = peak2NeutralMass - peak1NeutralMass;
                                                if (spacing1 > spacingRange.MaxValue || spacing1 < spacingRange.MinValue)
                                                {
                                                    pair.peaks[channelIndex + 3, j] = null;
                                                }
                                                if (spacing2 > spacingRange.MaxValue || spacing2 < spacingRange.MinValue)
                                                {
                                                    pair.peaks[channelIndex + 1, j] = null;
                                                    pair.peaks[channelIndex, j] = null;
                                                }
                                            }
                                        }
                                        else // peaks 2-4
                                        {
                                            peak2NeutralMass = Mass.MassFromMz(peak2.MZ, charge);
                                            peak3NeutralMass = Mass.MassFromMz(peak3.MZ, charge);
                                            peak4NeutralMass = Mass.MassFromMz(peak4.MZ, charge);
                                            spacing1 = peak4NeutralMass - peak3NeutralMass;
                                            spacing2 = peak3NeutralMass - peak2NeutralMass;
                                            if (spacing1 > spacingRange.MaxValue || spacing1 < spacingRange.MinValue)
                                            {
                                                pair.peaks[channelIndex + 3, j] = null;
                                            }
                                            if (spacing2 > spacingRange.MaxValue || spacing2 < spacingRange.MinValue)
                                            {
                                                pair.peaks[channelIndex + 2, j] = null;
                                                pair.peaks[channelIndex + 1, j] = null;
                                            }
                                        }
                                    }
                                    else
                                    {
                                        pair.peaks[channelIndex, j] = null;
                                        pair.peaks[channelIndex + 1, j] = null;
                                        pair.peaks[channelIndex + 2, j] = null;
                                        pair.peaks[channelIndex + 3, j] = null;
                                    }
                                }
                                else if (peakCount == 2)
                                {
                                    if (Form1.NOISEBANDCAP)
                                    {
                                        if (peak1 != null)
                                        {
                                            peak1NeutralMass = Mass.MassFromMz(peak1.MZ, charge);
                                            if (peak2 != null)
                                            {
                                                peak2NeutralMass = Mass.MassFromMz(peak2.MZ, charge);
                                                spacing1 = peak2NeutralMass - peak1NeutralMass;
                                            }
                                            else if (peak3 != null)
                                            {
                                                peak3NeutralMass = Mass.MassFromMz(peak3.MZ, charge);
                                                spacing1 = (peak3NeutralMass - peak1NeutralMass) / 2.0;
                                            }
                                            else
                                            {
                                                peak4NeutralMass = Mass.MassFromMz(peak4.MZ, charge);
                                                spacing1 = (peak4NeutralMass - peak1NeutralMass) / 3.0;
                                            }
                                        }
                                        else if (peak2 != null)
                                        {
                                            peak2NeutralMass = Mass.MassFromMz(peak2.MZ, charge);
                                            if (peak3 != null)
                                            {
                                                peak3NeutralMass = Mass.MassFromMz(peak3.MZ, charge);
                                                spacing1 = peak3NeutralMass - peak2NeutralMass;
                                            }
                                            else if (peak4 != null)
                                            {
                                                peak4NeutralMass = Mass.MassFromMz(peak4.MZ, charge);
                                                spacing1 = (peak4NeutralMass - peak2NeutralMass) / 2.0;
                                            }
                                            else
                                            {
                                                peak1NeutralMass = Mass.MassFromMz(peak1.MZ, charge);
                                                spacing1 = peak2NeutralMass - peak1NeutralMass;
                                            }
                                        }
                                        else if (peak3 != null)
                                        {
                                            peak3NeutralMass = Mass.MassFromMz(peak3.MZ, charge);
                                            if (peak4 != null)
                                            {
                                                peak4NeutralMass = Mass.MassFromMz(peak4.MZ, charge);
                                                spacing1 = peak4NeutralMass - peak3NeutralMass;
                                            }
                                            else if (peak1 != null)
                                            {
                                                peak1NeutralMass = Mass.MassFromMz(peak1.MZ, charge);
                                                spacing1 = (peak3NeutralMass - peak1NeutralMass) / 2.0;
                                            }
                                            else
                                            {
                                                peak2NeutralMass = Mass.MassFromMz(peak2.MZ, charge);
                                                spacing1 = peak3NeutralMass - peak2NeutralMass;
                                            }
                                        }
                                        else
                                        {
                                            peak4NeutralMass = Mass.MassFromMz(peak4.MZ, charge);
                                            if (peak1 != null)
                                            {
                                                peak1NeutralMass = Mass.MassFromMz(peak1.MZ, charge);
                                                spacing1 = (peak4NeutralMass - peak1NeutralMass) / 3.0;
                                            }
                                            else if (peak2 != null)
                                            {
                                                peak2NeutralMass = Mass.MassFromMz(peak2.MZ, charge);
                                                spacing1 = (peak4NeutralMass - peak2NeutralMass) / 2.0;
                                            }
                                            else
                                            {
                                                peak3NeutralMass = Mass.MassFromMz(peak3.MZ, charge);
                                                spacing1 = peak4NeutralMass - peak3NeutralMass;
                                            }
                                        }
                                        if (spacing1 > spacingRange.MaxValue || spacing1 < spacingRange.MinValue)
                                        {
                                            pair.peaks[channelIndex, j] = null;
                                            pair.peaks[channelIndex + 1, j] = null;
                                            pair.peaks[channelIndex + 2, j] = null;
                                            pair.peaks[channelIndex + 3, j] = null;
                                        }
                                    }
                                    else
                                    {
                                        pair.peaks[channelIndex, j] = null;
                                        pair.peaks[channelIndex + 1, j] = null;
                                        pair.peaks[channelIndex + 2, j] = null;
                                        pair.peaks[channelIndex + 3, j] = null;
                                    }
                                }
                                else // For 0 or 1 peaks, no further action needed
                                {
                                    pair.peaks[channelIndex, j] = null;
                                    pair.peaks[channelIndex + 1, j] = null;
                                    pair.peaks[channelIndex + 2, j] = null;
                                    pair.peaks[channelIndex + 3, j] = null;
                                }

                                //Check for # non-null peaks

                                bool stillQuantified = true;
                                for (int m = channelIndex; m < channelIndex + numIsotopologues; m++)
                                {
                                    if (pair.peaks[m, j] == null)
                                    {
                                        stillQuantified = false;
                                    }
                                }
                                if (stillQuantified)
                                {
                                    overallPeakCount++;
                                }
                            } // End isotopologue loop
                        } // End isotope loop
                    } // End cluster loop

                    quantified = false;
                    if (overallPeakCount > 0)
                    {
                        quantified = true;
                    }
                }
            }
            else
            {
                for (int i = 0; i < pairs.Count(); i++)
                {
                    pair = pairs[i];

                    for (int j = 0; j < numChannels; j += numIsotopologues)
                    {
                        for (int k = 0; k < numIsotopes; k++)
                        {
                            peak1 = pair.peaks[j, k];
                            peak2 = pair.peaks[j + 1, k];

                            //Only look at complete pairs
                            if (peak1 != null && peak2 != null)
                            {
                                peak1NeutralMass = Mass.MassFromMz(peak1.MZ, charge);
                                peak2NeutralMass = Mass.MassFromMz(peak2.MZ, charge);

                                spacing1 = peak2NeutralMass - peak1NeutralMass;
                                //theoSpacing = (double)lysineCount * lysineVersion;
                                //Spacing space = new Spacing(theoSpacing, spacing, charge, light.MZ, heavy.MZ);
                                //spacings.Add(space);
                                //Remove pairs whose spacing is more than 5 mDa away from the calculated spacing
                                if (spacing1 < spacingRange.MinValue || spacing1 > spacingRange.MaxValue)
                                {
                                    //Console.WriteLine("bad spacing");
                                    pair.peaks[j, k] = null;
                                    pair.peaks[j + 1, k] = null;
                                }
                            }
                            else if ((peak1 == null || peak2 == null) && Form1.NOISEBANDCAP)
                            {
                                // Do nothing
                            }
                            else
                            {
                                pair.peaks[j, k] = null;
                                pair.peaks[j + 1, k] = null;
                            }
                        }
                    }
                }
            }
        }

        public void applyNoise(ThermoRawFile rawFile)
        {
            PeptideSpectralMatch best = PSMs[rawFile][0];
            int charge = best.Charge;
            int scanNumber = best.ScanNumber;
            ThermoRawFileScan current;
            double injectionTime;
            try
            {
                if (Form1.MULTIINJECT)
                {
                    current = rawFile[scanNumber].GetNextMsnScan(2).GetNextMsnScan(1);
                }
                else
                {
                    current = rawFile[scanNumber].GetPrecursorRawScan();
                }
                injectionTime = rawFile.GetInjectionTime(current.ScanNum);

                double noise;
                List<double> noisePeaks = new List<double>();
                LabeledPeak closest;

                /*if (Form1.NEUCODE_SIXPLEX_ARG && arginineCount == 0)
                {
                    numChannels = 2;
                }*/

                //Calculate the average noise for all peaks found in the scan (includes closest peak for those not found in the scan)
                for (int c = 0; c < numChannels; c++)
                {
                    for (int j = 0; j < numIsotopes; j++)
                    {
                        if (peaks[c, j] != null)
                        {
                            noisePeaks.Add(peaks[c, j].Noise * injectionTime);
                        }
                        else
                        {
                            try
                            {
                                closest = (LabeledPeak)current.Spectrum.GetNearestPeak(Mass.MZFromMass(adjustedTheoMasses[c, j], charge));
                                noisePeaks.Add(closest.Noise * injectionTime);
                            }
                            catch (ArgumentOutOfRangeException)
                            {

                            }
                        }
                    }
                }
                noise = noisePeaks.Sum() / (double)noisePeaks.Count;

                //Apply noise to missing channels
                List<Pair> pairs;
                List<Pair> noNBCPairs;
                Pair pair;
                Pair noNBCPair;
                LabeledPeak peak;
                double coalescenceIntensity = 100000000.0;
                int channelIndex;

                if (coalescenceDetected)
                {
                    coalescedPeakIntensities.Sort();
                    coalescenceIntensity = coalescedPeakIntensities[0];
                }

                allHILACPairs.TryGetValue(rawFile, out pairs);
                noNBCHILACPairs.TryGetValue(rawFile, out noNBCPairs);
                for (int p = 0; p < pairs.Count(); p++)
                {
                    pair = pairs[p];
                    for (int j = 0; j < numIsotopes; j++)
                    {
                        for (int c = 0; c < numClusters; c++)
                        {
                            channelIndex = c * numIsotopologues;
                            //For complete pairs
                            if (pair.complete[c,j])
                            {
                                // If the pair falls above the peptide's coalescence threshold, set peaks to null
                                if (coalescenceDetected && pair.maxIntensity[c,j] > coalescenceIntensity)
                                {
                                    for (int m = channelIndex; m < channelIndex + numIsotopologues; m++)
                                    {
                                        pair.peaks[m, j] = null;
                                    }
                                }
                                // Otherwise, add pair to no NBC list
                                else
                                {
                                    bool found = false;
                                    int SN = pair.scanNumber;

                                    foreach (Pair noNBC in noNBCPairs)
                                    {
                                        if (noNBC.scanNumber == SN)
                                        {
                                            for (int m = channelIndex; m < channelIndex + numIsotopologues; m++)
                                            {
                                                noNBC.peaks[m, j] = pair.peaks[m, j];
                                            }
                                            found = true;
                                        }
                                    }

                                    if (!found)
                                    {
                                        noNBCPair = new Pair(this, rawFile, pair.scanNumber);
                                        for (int m = channelIndex; m < channelIndex + numIsotopologues; m++)
                                        {
                                            noNBCPair.peaks[m, j] = pair.peaks[m, j];
                                        }
                                        noNBCPairs.Add(noNBCPair);
                                    }
                                }
                            }
                            // For empty pairs
                            else if (pair.peakCount[c,j] == 0)
                            {
                                for (int m = 0; m < numChannels; m++)
                                {
                                    pair.peaks[m, j] = null;
                                }
                            }
                            // For incomplete pairs
                            else
                            {
                                // Apply noise to missing channels
                                if (coalescenceDetected && pair.maxIntensity[c,j] > coalescenceIntensity)
                                {
                                    if (pair.maxIntensity[c,j] > coalescenceIntensity)
                                    {
                                        for (int m = channelIndex; m < channelIndex + numIsotopologues; m++)
                                        {
                                            pair.peaks[m, j] = null;
                                        }
                                    }
                                }
                                else
                                {
                                    for (int m = channelIndex; m < channelIndex + numIsotopologues; m++)
                                    {
                                        if (pair.peaks[m, j] == null)
                                        {
                                            peak = new LabeledPeak();
                                            peak.MZ = 0;
                                            peak.dNL = noise;
                                            pair.peaks[m, j] = peak;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            catch (Exception)
            { }

        }

        //Calculates the systematic ppm error associated with the data set
        public void precursorPPMError(ThermoRawFile rawFile, List<PrecursorPPM> ppms)
        {
            PeptideSpectralMatch best = PSMs[rawFile].ElementAt(0);
            int charge = best.Charge;
            int scanNumber = best.ScanNumber;
            ThermoRawFileScan precursorScan;
            LabeledPeak peak;
            double precursorPPM;
            PrecursorPPM ppm;

            if (PSMs != null && PSMs.Count > 0)
            {
                //Match monoisotope peaks for systematic error calculation
                if (Form1.MULTIINJECT)
                {
                    precursorScan = rawFile[scanNumber].GetNextMsnScan(1).GetNextMsnScan(2);
                }
                else
                {
                    precursorScan = rawFile[scanNumber].GetPreviousMsnScan(1);
                }

                /*if (Form1.NEUCODE_SIXPLEX_ARG && arginineCount == 0)
                {
                    numChannels = 2;
                }
                else
                {
                    numChannels = NUMCHANNELS;
                }*/
                
                //Find all possible precursor peaks (monoisotope only)
                if (Form1.NEUCODE)
                {
                    double[] theo;
                    for (int c = 0; c < numClusters; c++)
                    {
                        theo = new double[numIsotopologues];
                        int channelIndex = c * numIsotopologues;
                        int count = 0;
                        for (int i = channelIndex; i < channelIndex + numIsotopologues; i++)
                        {
                            theo[count] = theoMasses[i, 0];
                            count++;
                        }
                        largestPeakPPM(theo, precursorScan, firstSearchRange, rawFile, ppms);
                    }
                }
                else
                {
                    if (Form1.NEUCODE_ARG_PROLINECONVERSION)
                    {
                        peak = largestPeak(theoMasses[0, 0], precursorScan, firstSearchRange, rawFile);
                        if (peak != null)
                        {
                            precursorPPM = Tolerance.GetError(Mass.MassFromMz(peak.MZ, charge), theoMasses[0, 0], ToleranceType.PPM);
                            ppm = new PrecursorPPM(charge, this.sequence, best.EValue, precursorPPM);
                            ppms.Add(ppm);
                        }
                    }
                    else
                    {
                        for (int i = 0; i < numChannels; i++)
                        {
                            peak = largestPeak(theoMasses[i, 0], precursorScan, firstSearchRange, rawFile);
                            if (peak != null)
                            {
                                precursorPPM = Tolerance.GetError(Mass.MassFromMz(peak.MZ, charge), theoMasses[i, 0], ToleranceType.PPM);
                                ppm = new PrecursorPPM(charge, this.sequence, best.EValue, precursorPPM);
                                ppms.Add(ppm);
                            }
                        }
                    }
                }
                
                
                /*if (Form1.NEUCODE_12PLEX)
                {
                    double[] theo;
                    for (int i = 0; i < numChannels; i += 4)
                    {
                        theo = new double[4];
                        if (Form1.FIRSTSEARCHDONE)
                        {
                            theo[0] = adjustedTheoMasses[i, 0];
                            theo[1] = adjustedTheoMasses[i + 1, 0];
                            theo[2] = adjustedTheoMasses[i + 2, 0];
                            theo[3] = adjustedTheoMasses[i + 3, 0];
                        }
                        else
                        {
                            theo[0] = theoMasses[i, 0];
                            theo[1] = theoMasses[i + 1, 0];
                            theo[2] = theoMasses[i + 2, 0];
                            theo[3] = theoMasses[i + 3, 0];
                        }
                        largestPeakPPM(theo, precursorScan, firstSearchRange, rawFile, ppms);
                    } 
                }*/
                /*else if (Form1.NEUCODE_DUPLEX_LYS1 || Form1.NEUCODE_DUPLEX_CARBAMYL)
                {
                    double[] theo;
                    for (int i = 0; i < numChannels; i += 2)
                    {
                        theo = new double[2];
                        if (Form1.FIRSTSEARCHDONE)
                        {
                            theo[0] = adjustedTheoMasses[i, 0];
                            theo[1] = adjustedTheoMasses[i + 1, 0];
                        }
                        else
                        {
                            theo[0] = theoMasses[i, 0];
                            theo[1] = theoMasses[i + 1, 0];
                        }
                        largestPeakPPM(theo, precursorScan, firstSearchRange, rawFile, ppms);
                    }
                }*/
                /*else if (Form1.NEUCODE_4PLEX_LIGHT || Form1.NEUCODE_4PLEX_HEAVY)
                {
                    double[] theo;
                    for (int i = 0; i < numChannels; i += 4)
                    {
                        theo = new double[4];
                        if (Form1.FIRSTSEARCHDONE)
                        {
                            theo[0] = adjustedTheoMasses[i, 0];
                            theo[1] = adjustedTheoMasses[i + 1, 0];
                            theo[2] = adjustedTheoMasses[i + 2, 0];
                            theo[3] = adjustedTheoMasses[i + 3, 0];
                        }
                        else
                        {
                            theo[0] = theoMasses[i, 0];
                            theo[1] = theoMasses[i + 1, 0];
                            theo[2] = theoMasses[i + 2, 0];
                            theo[3] = theoMasses[i + 3, 0];
                        }
                        largestPeakPPM(theo, precursorScan, firstSearchRange, rawFile, ppms);
                    }
                }*/
                /*else
                {
                    for (int i = 0; i < numChannels; i++)
                    {
                        peak = largestPeak(theoMasses[i, 0], precursorScan, firstSearchRange, rawFile);
                        if (peak != null)
                        {
                            precursorPPM = Tolerance.GetError(Mass.MassFromMz(peak.MZ, charge), theoMasses[i, 0], ToleranceType.PPM);
                            ppm = new PrecursorPPM(charge, this.sequence, best.EValue, precursorPPM);
                            ppms.Add(ppm);
                        }
                    }
                }*/             
            }
        }

        public double calculateMissingChannelFrequency(List<Pair> pairs)
        {
            double frequency = -1;
            int count = 0;
            int missingChannel = 0;

            if (pairs != null && pairs.Count > 0)
            {
                foreach (Pair pair in pairs)
                {
                    for (int c = 0; c < numChannels; c += numIsotopologues)
                    {
                        for (int j = 0; j < numIsotopes; j++)
                        {
                            if (numIsotopologues == 2)
                            {
                                if (pair.peaks[c, j] != null && pair.peaks[c + 1, j] != null)
                                {
                                    count++;
                                }
                                else if (pair.peaks[c, j] == null && pair.peaks[c + 1, j] == null)
                                {
                                    // Do not count pairs in which all peaks are null
                                }
                                else
                                {
                                    count++;
                                    missingChannel++;
                                }
                            }
                            if (numIsotopologues == 4)
                            {
                                if (pair.peaks[c, j] != null && pair.peaks[c+1,j] != null && pair.peaks[c+2,j] != null && pair.peaks[c+3,j] != null)
                                {
                                    count++;
                                }
                                else if (pair.peaks[c, j] == null && pair.peaks[c + 1, j] == null && pair.peaks[c + 2, j] == null && pair.peaks[c + 3, j] == null)
                                {
                                    // Do nothing when all peaks are null
                                }
                                else
                                {
                                    missingChannel++;
                                    count++;
                                }
                            }
                        }
                    }
                }
                frequency = (double)(missingChannel) / (double)(count);
            }
            return frequency;
        }

        public double[] calculateRatio(List<double> all, List<double> noNBC)
        {
            double[] statistics = new double[3];
            double allSum = 0;
            int allCount = all.Count();
            double allAverage;

            for (int i = 0; i < allCount; i++)
            {
                allSum += all[i];
            }
            allAverage = allSum / (double)allCount;

            double noNBCSum = 0;
            int noNBCCount = noNBC.Count();
            double noNBCAverage;

            for (int i = 0; i < noNBCCount; i++)
            {
                noNBCSum += noNBC[i];
            }
            noNBCAverage = noNBCSum / (double)noNBCCount;

            double allRSD = 0;
            double noNBCRSD = 0;

            for (int i = 0; i < allCount; i++)
            {
                allRSD += Math.Pow(all[i] - allAverage, 2.0);
            }
            allRSD = Math.Pow(allRSD / (double)allCount, 0.5) * 100.0;

            for (int i = 0; i < noNBCCount; i++)
            {
                noNBCRSD += Math.Pow(noNBC[i] - noNBCAverage, 2.0);
            }
            noNBCRSD = Math.Pow(noNBCRSD / (double)noNBCCount, 0.5) * 100.0;

            if (allRSD < noNBCRSD)
            {
                statistics[0] = allAverage;
                statistics[1] = allRSD;
                statistics[2] = allCount;
            }
            else
            {
                statistics[0] = noNBCAverage;
                statistics[1] = noNBCRSD;
                statistics[2] = noNBCCount;
            }
            return statistics;
        }

        public double[] calculateRatio(List<double> all)
        {
            double[] statistics = new double[3];
            double allSum = 0;
            int allCount = all.Count();
            double allAverage;

            for (int i = 0; i < allCount; i++)
            {
                allSum += all[i];
            }
            allAverage = allSum / (double)allCount;

            double allRSD = 0;

            for (int i = 0; i < allCount; i++)
            {
                allRSD += Math.Pow(all[i] - allAverage, 2.0);
            }
            allRSD = Math.Pow(allRSD / (double)allCount, 0.5) * 100.0;

            statistics[0] = allAverage;
            statistics[1] = allRSD;
            statistics[2] = allCount;

            return statistics;
        }

        public bool quantFilter(Pair pair, int isotope, int cluster, bool complete)
        {
            bool addIntensities = true;
            double[,] max;
            int channelIndex;

            /*if (Form1.NEUCODE_SIXPLEX_ARG && arginineCount == 0)
            {
                NUMCHANNELS = 2;
            }*/

            if (pair.peakCount[cluster, isotope] != numIsotopologues)
            {
                return false;
            }
            
            if (complete)
            {
                max = maxNoNBCIntensity;
                channelIndex = cluster * numIsotopologues;
                for (int i = channelIndex; i < channelIndex + numIsotopologues; i++)
                {
                    if (pair.peaks[i, isotope].dNL < INTENSITYCUTOFF * max[i, 0])
                    {
                        addIntensities = false;
                    }
                }
            }
            else
            {
                if (countNoNBCIsotopes[cluster] > 2)
                {
                    max = maxNoNBCIntensity;
                }
                else
                {
                    max = maxIntensity;
                }
                channelIndex = cluster * numIsotopologues;
                for (int i = channelIndex; i < channelIndex + numIsotopologues; i++)
                {
                    double mindNL = INTENSITYCUTOFF * max[i, 0];
                    if (pair.complete[cluster, isotope])
                    {
                        if (pair.peaks[i, isotope].dNL < mindNL)
                        {
                            addIntensities = false;
                        }
                    }
                    if (!pair.complete[cluster, isotope])
                    {
                        if (pair.peaks[i, isotope].MZ > 0 && pair.peaks[i, isotope].dNL < INTENSITYCUTOFF * max[i, 0])
                        {
                            addIntensities = false;
                        }
                    }
                }
            }
            return addIntensities;
        }

        public void quantify()
        {
            int[,] final = new int[numClusters, 2];
            quantifiedNBC = new bool[numClusters];
            finalQuantified = new int[numClusters];
            //totalIntensity = new double[NUMCHANNELS, NUMISOTOPES + 1];
            //noNBCTotalIntensity = new double[NUMCHANNELS, NUMISOTOPES + 1];
            //ratioList = new List<double>[(NUMCHANNELS / 2), 1];
            //noNBCRatioList = new List<double>[(NUMCHANNELS / 2), 1];
            //heavyToLightRatioSum = new double[NUMCHANNELS / 2, 1];
            
            if (noNBCHILACPairs != null && noNBCHILACPairs.Count > 0)
            {
                foreach (List<Pair> pairs in noNBCHILACPairs.Values)
                {
                    foreach (Pair pair in pairs)
                    {
                        for (int j = 0; j < numIsotopes; j++)
                        {
                            for (int c = 0; c < numClusters; c++)
                            {
                                if (Form1.QUANTFILTER)
                                {
                                    if (quantFilter(pair, j, c, true))
                                    {
                                        int channelIndex = c * numIsotopologues;
                                        for (int i = channelIndex; i < channelIndex + numIsotopologues; i++)
                                        {
                                            noNBCTotalIntensity[i, j] += pair.peaks[i, j].dNL;
                                            noNBCTotalIntensity[i, numIsotopes] += pair.peaks[i, j].dNL;

                                            /*if (noNBCRatioList[i / 2, 0] == null)
                                            {
                                                noNBCRatioList[i / 2, 0] = new List<double>();
                                            }
                                            noNBCRatioList[i / 2, 0].Add(pair.peaks[i+1,j].Intensity / pair.peaks[i,j].Intensity);*/
                                        }
                                        final[c, 1]++;
                                    }
                                }
                                else
                                {
                                    int channelIndex = c * numIsotopologues;
                                    bool noNullPeaks = true;
                                    for (int i = channelIndex; i < channelIndex + numIsotopologues; i++)
                                    {
                                        if (pair.peaks[i, j] == null)
                                        {
                                            noNullPeaks = false;
                                        }
                                    }

                                    if (noNullPeaks)
                                    {
                                        for (int i = channelIndex; i < channelIndex + numIsotopologues; i++)
                                        {
                                            noNBCTotalIntensity[i, j] += pair.peaks[i, j].dNL;
                                            noNBCTotalIntensity[i, numIsotopes] += pair.peaks[i, j].dNL;
                                        }
                                        final[c, 1]++;

                                        /*if (noNBCRatioList[i / 2, 0] == null)
                                        {
                                            noNBCRatioList[i / 2, 0] = new List<double>();
                                        }
                                        noNBCRatioList[i / 2, 0].Add(pair.peaks[i + 1, j].Intensity / pair.peaks[i, j].Intensity);*/
                                    }
                                }
                            } // End cluster loop
                        } // End isotope loop
                    } // End pair loop
                } // End pair list loop
            }
            if (allHILACPairs != null && allHILACPairs.Count > 0)
            {
                foreach (List<Pair> pairs in allHILACPairs.Values)
                {
                    foreach (Pair pair in pairs)
                    {
                        for (int j = 0; j < numIsotopes; j++)
                        {
                            for (int c = 0; c < numClusters; c++)
                            {
                                int channelIndex = c * numIsotopologues;
                                if (Form1.QUANTFILTER)
                                {
                                    //Use a peak's intensity if it is not noise-band capped and its intensity is greater than 1/2e of the maximum intensity

                                    if (quantFilter(pair, j, c, false))
                                    {
                                        for (int i = channelIndex; i < channelIndex + numIsotopologues; i++)
                                        {
                                            totalIntensity[i, j] += pair.peaks[i, j].dNL;
                                            totalIntensity[i, numIsotopes] += pair.peaks[i, j].dNL;

                                            /*if (ratioList[i / 2, 0] == null)
                                            {
                                                ratioList[i / 2, 0] = new List<double>();
                                            }
                                            ratioList[i / 2, 0].Add(pair.peaks[i + 1, j].Intensity / pair.peaks[i, j].Intensity);*/
                                        }
                                        final[c, 0]++;
                                    }
                                }
                                else
                                {
                                    bool noNullPeaks = true;
                                    int realPeaks = 0;
                                    for (int i = channelIndex; i < channelIndex + numIsotopologues; i++)
                                    {
                                        if (pair.peaks[i, j] == null)
                                        {
                                            noNullPeaks = false;
                                        }
                                        else
                                        {
                                            realPeaks++;
                                        }
                                    }

                                    if (noNullPeaks || (Form1.NEUCODE_ARG_PROLINECONVERSION && realPeaks > 0))
                                    {
                                        for (int i = channelIndex; i < channelIndex + numIsotopologues; i++)
                                        {
                                            if (pair.peaks[i, j] != null)
                                            {
                                                totalIntensity[i, j] += pair.peaks[i, j].dNL;
                                                totalIntensity[i, numIsotopes] += pair.peaks[i, j].dNL;
                                            }
                                        }
                                        final[c,0]++;

                                        /*if (noNBCRatioList[i / 2, 0] == null)
                                        {
                                            noNBCRatioList[i / 2, 0] = new List<double>();
                                        }
                                        noNBCRatioList[i / 2, 0].Add(pair.peaks[i + 1, j].Intensity / pair.peaks[i, j].Intensity);*/
                                    }
                                }
                            }
                        }
                    }
                }
            } 
            
            //heavyToLightRatioAverage = new double[NUMCHANNELS / 2, 1];
            //heavyToLightRatioMedian = new double[NUMCHANNELS / 2, 1];
            //heavyToLightRatioVariability = new double[NUMCHANNELS / 2, 1];
            //heavyToLightRatioCount = new double[NUMCHANNELS / 2, 1];

            double lightInt;
            double heavyInt;
            bool quantified = false;
            //int[,] ratioCount = new int[(NUMCHANNELS / 2), 1];
            //int[,] noNBCRatioCount = new int[(NUMCHANNELS / 2), 1];
            //double[,] ratioSum = new double[(NUMCHANNELS / 2), 1];
            //double[,] noNBCRatioSum = new double[(NUMCHANNELS / 2), 1];

            /*for (int i = 0; i < (NUMCHANNELS / 2); i++)
            {
                if (ratioList[i, 0] != null)
                {
                    ratioList[i, 0].Sort();
                    ratioCount[i, 0] = ratioList[i, 0].Count();

                    for (int j = 0; j < ratioCount[i, 0]; j++)
                    {
                        ratioSum[i, 0] += ratioList[i, 0][j];
                    }
                }

                if (noNBCRatioList[i, 0] != null)
                {
                    noNBCRatioList[i, 0].Sort();
                    noNBCRatioCount[i, 0] = noNBCRatioList[i, 0].Count();

                    for (int j = 0; j < noNBCRatioCount[i, 0]; j++)
                    {
                        noNBCRatioSum[i, 0] += noNBCRatioList[i, 0][j];
                    }
                }              
            }*/



            int minimumTotalPairs;
            int minimumNoNBCPairs;
            int minimumPostQFPairs;

            if (Form1.MULTIINJECT)
            {
                minimumTotalPairs = 1;
                minimumNoNBCPairs = 1;
                minimumPostQFPairs = 1;
            }
            else
            {
                minimumTotalPairs = 3;
                minimumNoNBCPairs = 3;
                minimumPostQFPairs = 3;
            }

            if (coalescenceDetected)
            {
                //Use only complete pairs for quantification
                for (int c = 0; c < numClusters; c++)
                {
                    int channelIndex = c * (numIsotopologues - 1);
                    for (int i = channelIndex + 1; i < channelIndex + numIsotopologues; i++)
                    {
                        int index1 = c * numIsotopologues;
                        int index2 = (c * numIsotopologues) + 1;
                        // Use only complete pairs for quantification
                        for (int n = channelIndex + 1; n < channelIndex + numIsotopologues; n++)
                        {
                            lightInt = noNBCTotalIntensity[index1, numIsotopes];
                            heavyInt = noNBCTotalIntensity[index2, numIsotopes];
                            finalQuantified[c] = final[c, 1];
                            if (lightInt > 0 && heavyInt > 0)
                            {
                                heavyToLightRatioSum[n - 1, 0] = heavyInt / lightInt;
                                quantifiedNBC[c] = false;
                            }
                            index2++;
                        }
                    }
                }
            }
            else
            {
                for (int c = 0; c < numClusters; c++)
                {
                    int channelIndex = c * (numIsotopologues - 1);

                    if (final[c, 1] >= minimumPostQFPairs) //&& ((double)countNoNBCIsotopes[c] / (double)countAllIsotopes[c]) > 0.25)
                    {
                        // Use only complete pairs for quantification
                        for (int i = channelIndex + 1; i < channelIndex + numIsotopologues; i++)
                        {
                            int index1 = c * numIsotopologues;
                            int index2 = (c * numIsotopologues) + 1;
                            // Use only complete pairs for quantification
                            for (int n = channelIndex + 1; n < channelIndex + numIsotopologues; n++)
                            {
                                lightInt = noNBCTotalIntensity[index1, numIsotopes];
                                heavyInt = noNBCTotalIntensity[index2, numIsotopes];
                                finalQuantified[c] = final[c, 1];
                                if (lightInt > 0 && heavyInt > 0)
                                {
                                    heavyToLightRatioSum[n - 1, 0] = heavyInt / lightInt;
                                    quantifiedNBC[c] = false;
                                }
                                index2++;
                            }
                        }
                    }
                    else if (final[c, 0] >= minimumPostQFPairs)
                    {
                        int index1 = c * numIsotopologues;
                        int index2 = (c * numIsotopologues) + 1;
                        // Use all pairs (complete & incomplete) for quantification
                        for (int i = channelIndex + 1; i < channelIndex + numIsotopologues; i++)
                        {
                            lightInt = totalIntensity[index1, numIsotopes];
                            heavyInt = totalIntensity[index2, numIsotopes];
                            finalQuantified[c] = final[c, 0];
                            if (lightInt > 0 && heavyInt > 0)
                            {
                                heavyToLightRatioSum[i - 1, 0] = heavyInt / lightInt;
                                quantifiedNBC[c] = true;
                            }
                            else if (Form1.NEUCODE_ARG_PROLINECONVERSION && lightInt > 0 && heavyInt == 0)
                            {
                                heavyToLightRatioSum[i - 1, 0] = 0;
                                quantifiedNBC[c] = true;
                            }
                            index2++;
                        }
                    }
                    else
                    {
                        // Not able to quantify
                        for (int i = channelIndex + 1; i < channelIndex + numIsotopologues; i++)
                        {
                            heavyToLightRatioSum[i - 1, 0] = double.NaN;
                        }
                    }
                }
            }
            /*if (countAllIsotopes >= minimumTotalPairs && countNoNBCIsotopes >= minimumNoNBCPairs && finalNBC >= minimumPostQFPairs  && !quantified && ((double) countNoNBCIsotopes / (double) countAllIsotopes) > 0.25)
            {
                // Use set of pairs with lowest % rsd for quantification
                for (int i = 1; i < numIsotopologues; i++)
                {
                    lightInt = noNBCTotalIntensity[0, numIsotopes];
                    heavyInt = noNBCTotalIntensity[i, numIsotopes];
                    finalQuantified = finalNBC;
                    if (lightInt > 0 && heavyInt > 0)
                    {
                        heavyToLightRatioSum[i-1, 0] = heavyInt / lightInt;
                        //heavyToLightRatioAverage[(i/2),0] = calculateRatio(ratioList[(i/2),0], noNBCRatioList[(i/2),0])[0];
                        //heavyToLightRatioVariability[(i / 2), 0] = calculateRatio(ratioList[(i / 2), 0], noNBCRatioList[(i / 2), 0])[1];
                        //heavyToLightRatioCount[(i / 2), 0] = calculateRatio(ratioList[(i / 2), 0], noNBCRatioList[(i / 2), 0])[2];
                        /*if (noNBCRatioCount[(i / 2), 0] % 2 != 0)
                        {
                            heavyToLightRatioMedian[(i / 2), 0] = noNBCRatioList[(i / 2), 0][noNBCRatioCount[(i / 2), 0] / 2];
                        }
                        else
                        {
                            heavyToLightRatioMedian[(i / 2), 0] = (noNBCRatioList[(i / 2), 0][noNBCRatioCount[(i / 2), 0] / 2] + noNBCRatioList[(i / 2), 0][(noNBCRatioCount[(i / 2), 0] / 2) - 1] )/ 2.0;
                        }

                        quantifiedNBC = false;
                        quantified = true;
                    }
                }
            }
            // Use all pairs for quantification
            if (countAllIsotopes >= minimumTotalPairs && !quantified)
            {
                for (int i = 1; i < numIsotopologues; i++)
                {
                    finalQuantified = final;
                    lightInt = totalIntensity[0, numIsotopes];
                    heavyInt = totalIntensity[i, numIsotopes];
                    if (lightInt > 0 && heavyInt > 0)
                    {
                        heavyToLightRatioSum[i - 1, 0] = heavyInt / lightInt;
                        quantifiedNBC = true;
                        //heavyToLightRatioAverage[(i / 2), 0] = calculateRatio(ratioList[(i / 2), 0])[0];
                        //heavyToLightRatioVariability[(i / 2), 0] = calculateRatio(ratioList[(i / 2), 0])[1];
                        //heavyToLightRatioCount[(i / 2), 0] = calculateRatio(ratioList[(i / 2), 0])[2];
                        
                        /*if (ratioCount[(i / 2), 0] % 2 != 0)
                        {
                            heavyToLightRatioMedian[(i / 2), 0] = ratioList[(i / 2), 0][ratioCount[(i / 2), 0] / 2];
                        }
                        else
                        {
                            heavyToLightRatioMedian[(i / 2), 0] = (ratioList[(i / 2), 0][ratioCount[(i / 2), 0] / 2] + ratioList[(i / 2), 0][(ratioCount[(i / 2), 0] / 2) - 1]) / 2.0;
                        }
                    }
                    else
                    {
                        heavyToLightRatioSum[i - 1, 0] = double.NaN;
                        //heavyToLightRatioAverage[(i / 2), 0] = double.NaN;
                        //heavyToLightRatioMedian[(i / 2), 0] = double.NaN;
                    }
                }
            }
            // Do not quantify
            if (countAllIsotopes < minimumTotalPairs && !quantified)
            {
                for (int i = 1; i < numIsotopologues; i++)
                {
                    heavyToLightRatioSum[i-1, 0] = double.NaN;
                    //heavyToLightRatioAverage[(i / 2), 0] = double.NaN;
                    //heavyToLightRatioMedian[(i / 2), 0] = double.NaN;
                }
            }*/
        }
    }
}
