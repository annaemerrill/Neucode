using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using CSMSL;
using CSMSL.Spectral;
using CSMSL.Chemistry;
using CSMSL.IO;
using CSMSL.Proteomics;
using CSMSL.IO.Thermo;

namespace Coon.NeuQuant
{
    class PeptideID
    {
        // Constants
        public static MassTolerance SILAC = new MassTolerance(MassToleranceType.PPM, 20.0); // Individual SILAC peak tolerance (for lower resolution MS1 scans)
        public static MassTolerance NEUCODE = new MassTolerance(MassToleranceType.PPM, 10.0); // Individual NeuCode peak tolerance (for higher resolution MS1 scans)
        public static double INTENSITYCUTOFF = 1.0 / (2.0 * Math.E); // Intensity threshold for quantitation filtering to eliminate low-level peaks
        public static double SIGNALTONOISE = Form1.MINIMUMSN;


        // Class members
        // Experimental design information
        public int numChannels; // # of total channels in the experiment (# clusters * # isotopologues)
        public int numIsotopes; // # of isotopes to include
        public int numClusters; // # of clusters
        public int numIsotopologues; // # of isotopologues within each cluster
        
        // Identification information
        public Dictionary<MSDataFile, List<PeptideSpectralMatch>> PSMs; // Organizes a peptide's PSMs based on raw file
        public PeptideSpectralMatch PSM; // A peptide's first PSM
        public Dictionary<MSDataFile, PeptideSpectralMatch> bestPSMs;
        public PeptideSpectralMatch bestPSM
        {
            get
            {
                PeptideSpectralMatch best = null;
                if (bestPSMs.Count > 0)
                {
                    foreach (PeptideSpectralMatch psm in bestPSMs.Values)
                    {
                        if (best == null || psm.EValue < best.EValue) best = psm;
                    }
                }
                return best;
            }
        }
        public int numLabels; // The number of quantitative labels carried by the peptide
        public string sequence; // A peptide's sequence
        public string scanNumbers
        {
            get
            {
                string scanNumbers = "";
                if (PSMs != null && PSMs.Count > 0)
                {
                    foreach (MSDataFile rawFile in PSMs.Keys) //For every raw file
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
        } // All the scan numbers that produced PSMs
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
        } // All the raw files that produced PSMs
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
        } // All the charge states that produced PSMs
        public double[,] theoMasses; // Theoretical masses of each peptide isotopologue 
        public double[,] adjustedTheoMasses
        {
            get
            {
                double ppmError = Form1.SYSTEMATICERROR;
                double adjustedMass;
                int channels;
                if (Form1.NEUCODE_SIXPLEX_ARG)
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
                                adjustedPrecursorMasses[i, j] = adjustedMass + (double)j * (Constants.CARBON13 - Constants.CARBON);
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
                            adjustedMass = theoMasses[i, 0];
                            adjustedPrecursorMasses[i, j] = adjustedMass + (double)j * (Constants.CARBON13 - Constants.CARBON);
                        }
                    }
                }
                return adjustedPrecursorMasses;
            }
        } // Theoretical masses of each peptide isotopologue after adjustment for systematic error
        public Range<double> rTRange; // Retention time window within which to look for MS1 scans
        public List<MSDataScan> fullScanList; // List of MS1 scans in which to look for peaks
        //public int bestPSMScanNumber
        //{
        //    get
        //    {
        //        int overallBestSN = 0;
        //        double lowestEValue = 100;
        //        if (PSMs != null && PSMs.Count > 0)
        //        {
        //            foreach (List<PeptideSpectralMatch> psms in PSMs.Values)
        //            {
        //                psms.OrderBy(PSM => PSM.EValue);
        //                if (psms[0].EValue < lowestEValue)
        //                {
        //                    lowestEValue = psms[0].EValue;
        //                    overallBestSN = psms[0].ScanNumber;
        //                }
        //            }
        //        }
        //        return overallBestSN;
        //    }
        //} // Finds the scan number with the best E-value
        
        // Quantification information
        public int peaksNeeded
        {
            get
            {
                int numPeaks;
                if (Form1.NOISEBANDCAP)
                {
                    if (!Form1.NEUCODE)
                    {
                        numPeaks = numChannels / 2;
                    }
                    else
                    {
                        numPeaks = numIsotopologues / 2;
                    }
                }
                else
                {
                    if (!Form1.NEUCODE)
                    {
                        numPeaks = numChannels;
                    }
                    else
                    {
                        numPeaks = numIsotopologues;
                    }
                }
                return numPeaks;
            }
        }
        public MassTolerance tolerance
        {
            get
            {
                MassTolerance tolerance;
                if (Form1.NEUCODE || Form1.SILAC_DUPLEX_LEUCN || Form1.SILAC_DUPLEX_LEUH)
                {
                    tolerance = NEUCODE;
                }
                else
                {
                    tolerance = SILAC;
                }
                return tolerance;
            }
        }
        public double[,] maxIntensity // Considering all pairs, finds the intensity and retention time at which each peptide isotopologue reaches its elution maximum
        {
            get
            {
                double[,] max = null;
                double injectionTime;
                MSDataFile rawFile;
                List<Pair> pairs;
                if (allPairs != null && allPairs.Count > 0)
                {
                    max = new double[numChannels, 2];
                    foreach (KeyValuePair<MSDataFile, List<Pair>> kvp in allPairs)
                    {
                        rawFile = kvp.Key;
                        pairs = kvp.Value;
                        foreach (Pair pair in pairs)
                        {
                            injectionTime = rawFile.GetInjectionTime(pair.scanNumber);
                            for (int i = 0; i < numChannels; i++)
                            {
                                for (int j = 0; j < numIsotopes; j++)
                                {
                                    if (pair.peaks[i,j] != null && pair.peaks[i, j].GetDenormalizedIntensity(injectionTime) > max[i, 0])
                                    {
                                        max[i, 0] = pair.peaks[i, j].GetDenormalizedIntensity(injectionTime);
                                        max[i, 1] = pair.rawFile[pair.scanNumber].RetentionTime; 
                                    }
                                }
                            }
                        }
                    }
                }
                return max;
            }
        }
        public double[,] maxCompleteIntensity // Considering complete pairs, finds the intensity and retention time at which each peptide isotopologue reaches its elution maximum
        {
            get
            {
                double[,] max = null;
                double injectionTime;
                MSDataFile rawFile;
                List<Pair> pairs;
                if (completePairs != null && completePairs.Count > 0)
                {
                    max = new double[numChannels, 2];
                    foreach (KeyValuePair<MSDataFile, List<Pair>> kvp in completePairs)
                    {
                        rawFile = kvp.Key;
                        pairs = kvp.Value;
                        foreach (Pair pair in pairs)
                        {
                            for (int i = 0; i < numChannels; i++)
                            {
                                for (int j = 0; j < numIsotopes; j++)
                                {
                                    if (pair.peaks[i, j] != null && pair.peaks[i, j].GetDenormalizedIntensity(pair.injectionTime) > max[i, 0])
                                    {
                                        max[i, 0] = pair.peaks[i, j].GetDenormalizedIntensity(pair.injectionTime);
                                        max[i, 1] = pair.rawFile[pair.scanNumber].RetentionTime;
                                    }
                                }
                            }
                        }
                    }
                }
                return max;
            }
        }
        public double maximumIntensity // Considering all isotopologues from all pairs, finds the peptide's maximum intensity
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
        public double maximumCompleteIntensity // Considering all isotopologues from complete pairs, finds the peptide's maximum intensity
        {
            get
            {
                double max = 0;
                if (maxCompleteIntensity != null)
                {
                    for (int c = 0; c < numChannels; c++)
                    {
                        if (maxCompleteIntensity[c, 0] > max)
                        {
                            max = maxCompleteIntensity[c, 0];
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
                List<Pair> all = new List<Pair>();
                double frequency = -1;
                if (allPairs != null && allPairs.Count > 0)
                {
                    foreach (List<Pair> pairs in allPairs.Values)
                    {
                        foreach (Pair pair in pairs)
                        {
                            all.Add(pair);
                        }
                    }
                    frequency = calculateMissingChannelFrequency(all);
                }
                return frequency;
            }
        } // Calculate's a peptide's missing channel frequency to be used in coalescence detection
        public double[,] completeTotalIntensity; // Keeps track of each channel's total intensity considering only complete pairs
        public double[,] totalIntensity; // Keeps track of each channel's total intensity considering all pairs
        public MassTolerance firstSearchMassRange; // The mass range used to search for peaks for systematic error determination
        public double[,] heavyToLightRatioSum; // Computes the quantitative ratio between channels
        public Dictionary<MSDataFile, List<Pair>> allPairs;
        public int countAllPairs
        {
            set { }
            get
            {
                int count = 0;
                if (allPairs != null && allPairs.Count > 0)
                {
                    foreach (MSDataFile rawFile in allPairs.Keys)
                    {
                        foreach (Pair pair in allPairs[rawFile])
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
                if (allPairs != null && allPairs.Count > 0)
                {
                    count = new int[numClusters];
                    foreach (MSDataFile rawFile in allPairs.Keys)
                    {
                        List<Pair> pairs = allPairs[rawFile];
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
                                    //if ((Form1.NEUCODE_ARG_PROLINECONVERSION && channelCount > 0) || channelCount == numIsotopologues)
                                    //{
                                    //    count[c]++;
                                    //}
                                }
                            }
                        }
                    }
                }
                return count;
            }
        }
        public Dictionary<MSDataFile, List<Pair>> completePairs;
        public int countCompletePairs
        {
            set { }
            get
            {
                int count = 0;
                if (completePairs != null && completePairs.Count > 0)
                {
                    foreach (MSDataFile rawFile in completePairs.Keys)
                    {
                        foreach (Pair pair in completePairs[rawFile])
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
        public int[] countCompleteIsotopes
        {
            set { }
            get
            {
                int[] count = null;
                int channelCount;
                int channelIndex;
                if (completePairs != null && completePairs.Count > 0)
                {
                    count = new int[numClusters];
                    foreach (MSDataFile rawFile in completePairs.Keys)
                    {
                        List<Pair> pairs = completePairs[rawFile];
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
        public bool[] quantifiedNoiseIncluded;
        public bool quantified;
        public int[] finalQuantified;
        public MassRange[] spacingMassRange
        {
            get
            {
                MassRange[] spacing;
                double theoSpacing1;
                double theoSpacing2;
                double theoSpacing3;
                double theoSpacing4;
                double theoSpacing5;
                if (numLabels > 0)
                {
                    if (Form1.NEUCODE)
                    {
                        spacing = new MassRange[numIsotopologues - 1];
                        if (numIsotopologues == 6)
                        {
                            theoSpacing1 = (double)numLabels * (2 * (Constants.DEUTERIUM - Constants.HYDROGEN) - 2 * (Constants.CARBON13 - Constants.CARBON));
                            theoSpacing2 = (double)numLabels * ((Constants.CARBON13 - Constants.CARBON) - (Constants.NITROGEN15 - Constants.NITROGEN));
                            theoSpacing3 = (double)numLabels * (4 * (Constants.DEUTERIUM - Constants.HYDROGEN) + (Constants.NITROGEN15 - Constants.NITROGEN) - 5 * (Constants.CARBON13 - Constants.CARBON));
                            theoSpacing4 = (double)numLabels * (4 * (Constants.CARBON13 - Constants.CARBON) - 2 * (Constants.DEUTERIUM - Constants.HYDROGEN) - 2 * (Constants.NITROGEN15 - Constants.NITROGEN));
                            theoSpacing5 = (double)numLabels * (4 * (Constants.DEUTERIUM - Constants.HYDROGEN) - 4 * (Constants.CARBON13 - Constants.CARBON));
                            spacing[0] = new MassRange(theoSpacing1, new MassTolerance(MassToleranceType.DA, 0.5 * theoSpacing1));
                            spacing[1] = new MassRange(theoSpacing2, new MassTolerance(MassToleranceType.DA, 0.5 * theoSpacing2));
                            spacing[2] = new MassRange(theoSpacing3, new MassTolerance(MassToleranceType.DA, 0.5 * theoSpacing3));
                            spacing[3] = new MassRange(theoSpacing4, new MassTolerance(MassToleranceType.DA, 0.5 * theoSpacing4));
                            spacing[4] = new MassRange(theoSpacing5, new MassTolerance(MassToleranceType.DA, 0.5 * theoSpacing5));
                        }
                        else if (numIsotopologues == 4)
                        {
                            if (Form1.NEUCODE_FOURPLEX_LYS8_12MDA)
                            {
                                theoSpacing1 = (double)numLabels * (2 * (Constants.DEUTERIUM - Constants.HYDROGEN) - (Constants.NITROGEN15 - Constants.NITROGEN) - (Constants.CARBON13 - Constants.CARBON));
                                theoSpacing2 = theoSpacing1;
                                theoSpacing3 = (double)numLabels * (4 * (Constants.DEUTERIUM - Constants.HYDROGEN) - 4 * (Constants.CARBON13 - Constants.CARBON));
                                spacing[0] = new MassRange(theoSpacing1, new MassTolerance(MassToleranceType.DA, 0.5 * theoSpacing1));
                                spacing[1] = spacing[0];
                                spacing[2] = new MassRange(theoSpacing3, new MassTolerance(MassToleranceType.DA, 0.5 * theoSpacing3));
                            }
                            else
                            {
                                theoSpacing1 = (double)numLabels * (2 * (Constants.CARBON13 - Constants.CARBON) - 2 * (Constants.NITROGEN15 - Constants.NITROGEN));
                                theoSpacing2 = theoSpacing1;
                                theoSpacing3 = theoSpacing1;
                                spacing[0] = new MassRange(theoSpacing1, new MassTolerance(MassToleranceType.DA, 0.5 * theoSpacing1));
                                spacing[1] = spacing[0];
                                spacing[2] = spacing[1];
                            }
                        }
                        else if (numIsotopologues == 3)
                        {
                            theoSpacing1 = (double)numLabels * (6 * (Constants.DEUTERIUM - Constants.HYDROGEN) - 6 * (Constants.CARBON13 - Constants.CARBON));
                            theoSpacing2 = (double)numLabels * (2 * (Constants.DEUTERIUM - Constants.HYDROGEN) - 2 * (Constants.NITROGEN15 - Constants.NITROGEN));
                            spacing[0] = new MassRange(theoSpacing1, new MassTolerance(MassToleranceType.DA, 0.5 * theoSpacing1));
                            spacing[1] = new MassRange(theoSpacing2, new MassTolerance(MassToleranceType.DA, 0.5 * theoSpacing2));
                        }
                        else if (Form1.NEUCODE_DUPLEX_LYS1_6MDA || Form1.NEUCODE_DUPLEX_CARBAMYL)
                        {
                            theoSpacing1 = (double)numLabels * ((Constants.CARBON13 - Constants.CARBON) - (Constants.NITROGEN15 - Constants.NITROGEN));
                            spacing[0] = new MassRange(theoSpacing1, new MassTolerance(MassToleranceType.DA, 0.5 * theoSpacing1));
                        }
                        else if (Form1.NEUCODE_DUPLEX_LEU7_18MDA)
                        {
                            if (conversionFactor > 0)
                            {
                                theoSpacing1 = ((double)numLabels * (7 * (Constants.DEUTERIUM - Constants.HYDROGEN) - (6 * (Constants.CARBON13 - Constants.CARBON) + (1 * (Constants.NITROGEN15 - Constants.NITROGEN))))) - conversionFactor * (Constants.NITROGEN15 - Constants.NITROGEN);
                            }
                            else
                            {
                                theoSpacing1 = (double)numLabels * (7 * (Constants.DEUTERIUM - Constants.HYDROGEN) - (6 * (Constants.CARBON13 - Constants.CARBON) + (1 * (Constants.NITROGEN15 - Constants.NITROGEN))));
                            }
                            spacing[0] = new MassRange(theoSpacing1, new MassTolerance(MassToleranceType.DA, 0.5 * theoSpacing1));
                        }
                        else
                        {
                            theoSpacing1 = (double)numLabels * (8 * (Constants.DEUTERIUM - Constants.HYDROGEN) - (6 * (Constants.CARBON13 - Constants.CARBON) + (2 * (Constants.NITROGEN15 - Constants.NITROGEN))));
                            spacing[0] = new MassRange(theoSpacing1, new MassTolerance(MassToleranceType.DA, 0.010 * numLabels));
                        }

                    }
                    else
                    {
                        spacing = new MassRange[numIsotopologues];
                        if (Form1.SILAC_DUPLEX_LYSH)
                        {
                            theoSpacing1 = numLabels * (8 * (Constants.DEUTERIUM - Constants.HYDROGEN));
                        }
                        else if (Form1.SILAC_DUPLEX_LEUCN)
                        {
                            if (conversionFactor > 0)
                            {
                                theoSpacing1 = numLabels * ((6 * (Constants.CARBON13 - Constants.CARBON)) + (Constants.NITROGEN15 - Constants.NITROGEN)) - conversionFactor * (Constants.NITROGEN15 - Constants.NITROGEN);
                            }
                            else
                            {
                                theoSpacing1 = numLabels * ((6 * (Constants.CARBON13 - Constants.CARBON)) + (Constants.NITROGEN15 - Constants.NITROGEN));
                            }
                        }
                        else if (Form1.SILAC_DUPLEX_LEUH)
                        {
                            theoSpacing1 = numLabels * (7 * (Constants.DEUTERIUM - Constants.HYDROGEN));
                        }
                        else
                        {
                            theoSpacing1 = numLabels * ((6 * (Constants.CARBON13 - Constants.CARBON)) + (2 * (Constants.NITROGEN15 - Constants.NITROGEN)));
                        }
                        spacing[0] = new MassRange(theoSpacing1, new MassTolerance(MassToleranceType.DA, 0.010 * numLabels));
                    }
                    return spacing;
                }
                else
                {
                    return null;
                }
            }
        }
        public int conversionFactor { get; set; }
        public bool theoreticallyResolvable
        {
            get
            {
                if (!Form1.NEUCODE)
                {
                    return true;
                }
                else
                {
                    int charge = bestPSM.Charge;
                    double FWTM = 1.82262 * ((Mass.MzFromMass(theoMasses[0,0], charge) / (480000.0 * Math.Sqrt(400 / Mass.MzFromMass(theoMasses[0,0], charge)))));
                    double separation = (numLabels * spacingMassRange[0].Mean) / (double)charge;

                    if (separation > FWTM)
                    {
                        return true;
                    }
                    return false;
                }
            }
        }

        // Coalescence information
        public bool coalescenceDetected;
        public List<double> coalescedPeakIntensities;

        // Creates a PeptideID based on a database search peptide identification
        public PeptideID(int scanNumber, int charge, double eValue, string sequenceOriginal, MSDataFile rawFile, string mods)
        {           
            // Local variables
            List<PeptideSpectralMatch> psms;
            Peptide peptide;
            Peptide check1;
            Peptide check2;
            Peptide check3;
            Peptide check4;
            Peptide check5;
            Peptide check6;
            Peptide check7;
            Peptide check8;
            Peptide check9;
            Peptide check10;
            Peptide check11;
            Peptide check12;

            // Get rid of variable label incorporations
            string sequenceFixed = "";
            for (int i = 0; i < sequenceOriginal.Length; i++)
            {
                if (sequenceOriginal[i].Equals('l')) sequenceFixed += 'L';
                else sequenceFixed += sequenceOriginal[i];
            }
            
            // Initialize all peptide properties
            sequence = sequenceFixed;
            numChannels = Form1.NUMCHANNELS;
            numIsotopes = Form1.NUMISOTOPES;
            numIsotopologues = Form1.NUMISOTOPOLOGUES;
            numClusters = Form1.NUMCLUSTERS;

            // Deal with variable mods that do not affect quantification (e.g., oxidation, phosphorylation)
            string sequenceNoMods = "";
            List<int> oxidationPositions = new List<int>();
            List<int> phosphorylationPositions = new List<int>();
            List<int> tyrosineNHSPositions = new List<int>();
            for (int i = 0; i < sequence.Length; i++)
            {
                if (sequence[i].Equals('m'))
                {
                    sequenceNoMods += 'M';
                    oxidationPositions.Add(i);
                }
                else if (sequence[i].Equals('s'))
                {
                    sequenceNoMods += 'S';
                    phosphorylationPositions.Add(i);
                }
                else if (sequence[i].Equals('t'))
                {
                    sequenceNoMods += 'T';
                    phosphorylationPositions.Add(i);
                }
                else if (sequence[i].Equals('y'))
                {
                    sequenceNoMods += 'Y';
                    if (Form1.NEUCODE_DUPLEX_CARBAMYL || Form1.NEUCODE_4PLEX_LIGHT || Form1.NEUCODE_4PLEX_MEDIUM || Form1.NEUCODE_4PLEX_HEAVY || Form1.NEUCODE_12PLEX)
                    {
                        tyrosineNHSPositions.Add(i);
                    }
                    else
                    {
                        phosphorylationPositions.Add(i);
                    }
                }
                else if (sequence[i].Equals('k'))
                {
                    sequenceNoMods += 'K';
                }
                else
                {
                    sequenceNoMods += sequence[i];
                }
            }

            peptide = new Peptide(sequenceNoMods);
            peptide.SetModification(NamedChemicalFormula.Carbamidomethyl, 'C');
            if (oxidationPositions.Count > 0)
            {
                foreach (int position in oxidationPositions)
                {
                    peptide.SetModification(NamedChemicalFormula.Oxidation, position + 1);
                }
            }
            if (phosphorylationPositions.Count > 0)
            {
                foreach (int position in phosphorylationPositions)
                {
                    peptide.SetModification(NamedChemicalFormula.Phosphorylation, position);
                }
            }

            // Set number of labels and search range for PPM correction
            if (Form1.NEUCODE)
            {
                if (Form1.NEUCODE_DUPLEX_CARBAMYL || Form1.NEUCODE_4PLEX_LIGHT || Form1.NEUCODE_4PLEX_MEDIUM || Form1.NEUCODE_4PLEX_HEAVY || Form1.NEUCODE_12PLEX)
                {
                    numLabels = countResidues('K', peptide.Sequence) + 1 + tyrosineNHSPositions.Count;
                    firstSearchMassRange = new MassTolerance(MassToleranceType.PPM, 50.0);
                }
                else if (Form1.NEUCODE_DUPLEX_LEU7_18MDA)
                {
                    numLabels = countResidues('L', peptide.Sequence);
                    firstSearchMassRange = new MassTolerance(MassToleranceType.PPM, 50.0);
                }
                else
                {
                    numLabels = countResidues('K', peptide.Sequence);
                    firstSearchMassRange = new MassTolerance(MassToleranceType.PPM, 50.0);
                }
            }
            else
            {
                if (Form1.SILAC_DUPLEX_LEUCN || Form1.SILAC_DUPLEX_LEUH)
                {
                    numLabels = countResidues('L', peptide.Sequence);
                    firstSearchMassRange = new MassTolerance(MassToleranceType.PPM, 50.0);
                }
                else
                {
                    numLabels = countResidues('K', peptide.Sequence);
                    firstSearchMassRange = new MassTolerance(MassToleranceType.PPM, 100.0);
                }
            }

            theoMasses = new double[numChannels, 1];

            // Initialize all identification data members
            PSMs = new Dictionary<MSDataFile, List<PeptideSpectralMatch>>();
            PSMs.Add(rawFile, new List<PeptideSpectralMatch>());
            PSM = new PeptideSpectralMatch(scanNumber, charge, eValue);
            PSMs.TryGetValue(rawFile, out psms);
            psms.Add(PSM);
            bestPSMs = new Dictionary<MSDataFile, PeptideSpectralMatch>();
            bestPSMs.Add(rawFile, PSM);
            
            // Initialize all quantitation data members
            completeTotalIntensity = new double[numChannels, numIsotopes + 1];
            totalIntensity = new double[numChannels, numIsotopes + 1];
            if (Form1.NEUCODE)
            {
                heavyToLightRatioSum = new double[(numIsotopologues - 1) * numClusters, 1];
            }
            else
            {
                heavyToLightRatioSum = new double[(numChannels - 1), 1];
            }
            allPairs = new Dictionary<MSDataFile, List<Pair>>();
            allPairs.Add(rawFile, new List<Pair>());
            completePairs = new Dictionary<MSDataFile, List<Pair>>();
            completePairs.Add(rawFile, new List<Pair>()); 
            
            // Set theoretical masses of each channel
            if (numLabels > 0)
            {             
                if (Form1.NEUCODE_DUPLEX_LYS8_36MDA || Form1.NEUCODE_DUPLEX_LYS1_6MDA || Form1.NEUCODE_DUPLEX_LEU7_18MDA)
                {
                    check1 = new Peptide(peptide);
                    check2 = new Peptide(peptide);
                    if (Form1.NEUCODE_DUPLEX_LYS8_36MDA)
                    {
                        check1.SetModification(NamedChemicalFormula.GetModification("Lys +8 13C6 15N2"), ModificationSites.K);
                        check2.SetModification(NamedChemicalFormula.GetModification("Lys +8 2H8"), ModificationSites.K);
                    }
                    else if (Form1.NEUCODE_DUPLEX_LYS1_6MDA)
                    {
                        check1.SetModification(NamedChemicalFormula.GetModification("Lys +1 15N"), ModificationSites.K);
                        check2.SetModification(NamedChemicalFormula.GetModification("Lys +1 13C"), ModificationSites.K);
                    }
                    theoMasses[0, 0] = check1.Mass.Monoisotopic;
                    theoMasses[1, 0] = check2.Mass.Monoisotopic;
                }

                else if (Form1.NEUCODE_TRIPLEX_LYS8_18MDA)
                {
                    check1 = new Peptide(peptide);
                    check2 = new Peptide(peptide);
                    check3 = new Peptide(peptide);
                    check1.SetModification(NamedChemicalFormula.GetModification("Lys +8 13C6 15N2"), ModificationSites.K);
                    check2.SetModification(NamedChemicalFormula.GetModification("Lys +8 2H6 15N2"), ModificationSites.K);
                    check3.SetModification(NamedChemicalFormula.GetModification("Lys +8 2H8"), ModificationSites.K);
                    theoMasses[0, 0] = check1.Mass.Monoisotopic;
                    theoMasses[1, 0] = check2.Mass.Monoisotopic;
                    theoMasses[2, 0] = check3.Mass.Monoisotopic;
                }

                else if (Form1.NEUCODE_FOURPLEX_LYS8_12MDA)
                {
                    check1 = new Peptide(peptide);
                    check2 = new Peptide(peptide);
                    check3 = new Peptide(peptide);
                    check4 = new Peptide(peptide);
                    check1.SetModification(NamedChemicalFormula.GetModification("Lys +8 13C6 15N2"), ModificationSites.K);
                    check2.SetModification(NamedChemicalFormula.GetModification("Lys +8 13C5 2H2 15N1"), ModificationSites.K);
                    check3.SetModification(NamedChemicalFormula.GetModification("Lys +8 13C4 2H4"), ModificationSites.K);
                    check4.SetModification(NamedChemicalFormula.GetModification("Lys +8 2H8"), ModificationSites.K);
                    theoMasses[0, 0] = check1.Mass.Monoisotopic;
                    theoMasses[1, 0] = check2.Mass.Monoisotopic;
                    theoMasses[2, 0] = check3.Mass.Monoisotopic;
                    theoMasses[3, 0] = check4.Mass.Monoisotopic;
                }

                else if (Form1.NEUCODE_SIXPLEX_LYS8_6MDA)
                {
                    check1 = new Peptide(peptide);
                    check2 = new Peptide(peptide);
                    check3 = new Peptide(peptide);
                    check4 = new Peptide(peptide);
                    check5 = new Peptide(peptide);
                    check6 = new Peptide(peptide);
                    check1.SetModification(NamedChemicalFormula.GetModification("Lys +8 13C6 15N2"), ModificationSites.K);
                    check2.SetModification(NamedChemicalFormula.GetModification("Lys +8 13C4 2H2 15N2"), ModificationSites.K);
                    check3.SetModification(NamedChemicalFormula.GetModification("Lys +8 13C5 2H2 15N1"), ModificationSites.K);
                    check4.SetModification(NamedChemicalFormula.GetModification("Lys +8 2H6 15N2"), ModificationSites.K);
                    check5.SetModification(NamedChemicalFormula.GetModification("Lys +8 13C4 2H4"), ModificationSites.K);
                    check6.SetModification(NamedChemicalFormula.GetModification("Lys +8 2H8"), ModificationSites.K);
                    theoMasses[0, 0] = check1.Mass.Monoisotopic;
                    theoMasses[1, 0] = check2.Mass.Monoisotopic;
                    theoMasses[2, 0] = check3.Mass.Monoisotopic;
                    theoMasses[3, 0] = check4.Mass.Monoisotopic;
                    theoMasses[4, 0] = check5.Mass.Monoisotopic;
                    theoMasses[5, 0] = check6.Mass.Monoisotopic;
                }

                else if (Form1.NEUCODE_DUPLEX_CARBAMYL)
                {
                    check1 = new Peptide(peptide);
                    check2 = new Peptide(peptide);

                    check1.SetModification(NamedChemicalFormula.GetModification("Carbamyl L"), ModificationSites.K | ModificationSites.NPep);
                    check2.SetModification(NamedChemicalFormula.GetModification("Carbamyl H"), ModificationSites.K | ModificationSites.NPep);
                    if (tyrosineNHSPositions.Count > 0)
                    {
                        foreach (int position in tyrosineNHSPositions)
                        {
                            check1.SetModification(NamedChemicalFormula.GetModification("Carbamyl L"), position);
                            check2.SetModification(NamedChemicalFormula.GetModification("Carbamyl H"), position);
                        }
                    }
                    theoMasses[0, 0] = check1.Mass.Monoisotopic;
                    theoMasses[1, 0] = check2.Mass.Monoisotopic;
                }

                else if (Form1.NEUCODE_4PLEX_LIGHT || Form1.NEUCODE_4PLEX_MEDIUM || Form1.NEUCODE_4PLEX_HEAVY)
                {
                    check1 = new Peptide(peptide);
                    check2 = new Peptide(peptide);
                    check3 = new Peptide(peptide);
                    check4 = new Peptide(peptide);
                    string mod1;
                    string mod2;
                    string mod3;
                    string mod4;
                    if (Form1.NEUCODE_4PLEX_LIGHT)
                    {
                        mod1 = "4plex L1";
                        mod2 = "4plex L2";
                        mod3 = "4plex L3";
                        mod4 = "4plex L4";
                    }
                    else if (Form1.NEUCODE_4PLEX_MEDIUM)
                    {
                        mod1 = "4plex M1";
                        mod2 = "4plex M2";
                        mod3 = "4plex M3";
                        mod4 = "4plex M4";
                    }
                    else
                    {
                        mod1 = "4plex H1";
                        mod2 = "4plex H2";
                        mod3 = "4plex H3";
                        mod4 = "4plex H4";
                    }
                    check1.SetModification(NamedChemicalFormula.GetModification(mod1), ModificationSites.K | ModificationSites.NPep);
                    check2.SetModification(NamedChemicalFormula.GetModification(mod2), ModificationSites.K | ModificationSites.NPep);
                    check3.SetModification(NamedChemicalFormula.GetModification(mod3), ModificationSites.K | ModificationSites.NPep);
                    check4.SetModification(NamedChemicalFormula.GetModification(mod4), ModificationSites.K | ModificationSites.NPep);
                    if (tyrosineNHSPositions.Count > 0)
                    {
                        foreach (int position in tyrosineNHSPositions)
                        {
                            check1.SetModification(NamedChemicalFormula.GetModification(mod1), position);
                            check2.SetModification(NamedChemicalFormula.GetModification(mod2), position);
                            check3.SetModification(NamedChemicalFormula.GetModification(mod3), position);
                            check4.SetModification(NamedChemicalFormula.GetModification(mod4), position);
                        }
                    }
                    theoMasses[0, 0] = check1.Mass.Monoisotopic;
                    theoMasses[1, 0] = check2.Mass.Monoisotopic;
                    theoMasses[2, 0] = check3.Mass.Monoisotopic;
                    theoMasses[3, 0] = check4.Mass.Monoisotopic;
                }

                else if (Form1.NEUCODE_12PLEX)
                {
                    check1 = new Peptide(peptide);
                    check2 = new Peptide(peptide);
                    check3 = new Peptide(peptide);
                    check4 = new Peptide(peptide);
                    check5 = new Peptide(peptide);
                    check6 = new Peptide(peptide);
                    check7 = new Peptide(peptide);
                    check8 = new Peptide(peptide);
                    check9 = new Peptide(peptide);
                    check10 = new Peptide(peptide);
                    check11 = new Peptide(peptide);
                    check12 = new Peptide(peptide);

                    check1.SetModification(NamedChemicalFormula.GetModification("4plex L1"), ModificationSites.K | ModificationSites.NPep);
                    check2.SetModification(NamedChemicalFormula.GetModification("4plex L2"), ModificationSites.K | ModificationSites.NPep);
                    check3.SetModification(NamedChemicalFormula.GetModification("4plex L3"), ModificationSites.K | ModificationSites.NPep);
                    check4.SetModification(NamedChemicalFormula.GetModification("4plex L4"), ModificationSites.K | ModificationSites.NPep);
                    check5.SetModification(NamedChemicalFormula.GetModification("4plex M1"), ModificationSites.K | ModificationSites.NPep);
                    check6.SetModification(NamedChemicalFormula.GetModification("4plex M2"), ModificationSites.K | ModificationSites.NPep);
                    check7.SetModification(NamedChemicalFormula.GetModification("4plex M3"), ModificationSites.K | ModificationSites.NPep);
                    check8.SetModification(NamedChemicalFormula.GetModification("4plex M4"), ModificationSites.K | ModificationSites.NPep);
                    check9.SetModification(NamedChemicalFormula.GetModification("4plex H1"), ModificationSites.K | ModificationSites.NPep);
                    check10.SetModification(NamedChemicalFormula.GetModification("4plex H2"), ModificationSites.K | ModificationSites.NPep);
                    check11.SetModification(NamedChemicalFormula.GetModification("4plex H3"), ModificationSites.K | ModificationSites.NPep);
                    check12.SetModification(NamedChemicalFormula.GetModification("4plex H4"), ModificationSites.K | ModificationSites.NPep);
                    if (tyrosineNHSPositions.Count > 0)
                    {
                        foreach (int position in tyrosineNHSPositions)
                        {
                            check1.SetModification(NamedChemicalFormula.GetModification("4plex L1"), position);
                            check2.SetModification(NamedChemicalFormula.GetModification("4plex L2"), position);
                            check3.SetModification(NamedChemicalFormula.GetModification("4plex L3"), position);
                            check4.SetModification(NamedChemicalFormula.GetModification("4plex L4"), position);
                            check5.SetModification(NamedChemicalFormula.GetModification("4plex M1"), position);
                            check6.SetModification(NamedChemicalFormula.GetModification("4plex M2"), position);
                            check7.SetModification(NamedChemicalFormula.GetModification("4plex M3"), position);
                            check8.SetModification(NamedChemicalFormula.GetModification("4plex M4"), position);
                            check9.SetModification(NamedChemicalFormula.GetModification("4plex H1"), position);
                            check10.SetModification(NamedChemicalFormula.GetModification("4plex H2"), position);
                            check11.SetModification(NamedChemicalFormula.GetModification("4plex H3"), position);
                            check12.SetModification(NamedChemicalFormula.GetModification("4plex H4"), position);
                        }
                    }

                    theoMasses[0, 0] = check1.Mass.Monoisotopic;
                    theoMasses[1, 0] = check2.Mass.Monoisotopic;
                    theoMasses[2, 0] = check3.Mass.Monoisotopic;
                    theoMasses[3, 0] = check4.Mass.Monoisotopic;
                    theoMasses[4, 0] = check5.Mass.Monoisotopic;
                    theoMasses[5, 0] = check6.Mass.Monoisotopic;
                    theoMasses[6, 0] = check7.Mass.Monoisotopic;
                    theoMasses[7, 0] = check8.Mass.Monoisotopic;
                    theoMasses[8, 0] = check9.Mass.Monoisotopic;
                    theoMasses[9, 0] = check10.Mass.Monoisotopic;
                    theoMasses[10, 0] = check11.Mass.Monoisotopic;
                    theoMasses[11, 0] = check12.Mass.Monoisotopic;
                }

                else if (Form1.NEUCODE_SIXPLEX_MTRAQ)
                {
                    check1 = new Peptide(peptide);
                    check2 = new Peptide(peptide);
                    check3 = new Peptide(peptide);
                    check4 = new Peptide(peptide);
                    check5 = new Peptide(peptide);
                    check6 = new Peptide(peptide);

                    check1.SetModification(NamedChemicalFormula.GetModification("mTRAQ L Lys +8 13C6 15N2"), ModificationSites.K | ModificationSites.NPep);
                    check2.SetModification(NamedChemicalFormula.GetModification("mTRAQ L Lys +8 2H8"), ModificationSites.K | ModificationSites.NPep);
                    check3.SetModification(NamedChemicalFormula.GetModification("mTRAQ M Lys +8 13C6 15N2"), ModificationSites.K | ModificationSites.NPep);
                    check4.SetModification(NamedChemicalFormula.GetModification("mTRAQ M Lys +8 2H8"), ModificationSites.K | ModificationSites.NPep);
                    check5.SetModification(NamedChemicalFormula.GetModification("mTRAQ H Lys +8 13C6 15N2"), ModificationSites.K | ModificationSites.NPep);
                    check6.SetModification(NamedChemicalFormula.GetModification("mTRAQ H Lys +8 2H8"), ModificationSites.K | ModificationSites.NPep);

                    theoMasses[0, 0] = check1.Mass.Monoisotopic;
                    theoMasses[1, 0] = check2.Mass.Monoisotopic;
                    theoMasses[2, 0] = check3.Mass.Monoisotopic;
                    theoMasses[3, 0] = check4.Mass.Monoisotopic;
                    theoMasses[4, 0] = check5.Mass.Monoisotopic;
                    theoMasses[5, 0] = check6.Mass.Monoisotopic;
                }

                else if (Form1.NEUCODE_SIXPLEX_ARG || Form1.NEUCODE_SIXPLEX_LEU)
                {
                    check1 = new Peptide(peptide);
                    check2 = new Peptide(peptide);
                    check3 = new Peptide(peptide);
                    check4 = new Peptide(peptide);
                    check5 = new Peptide(peptide);
                    check6 = new Peptide(peptide);

                    check1.SetModification(NamedChemicalFormula.GetModification("Lys +8 13C6 15N2"), ModificationSites.K);
                    check2.SetModification(NamedChemicalFormula.GetModification("Lys +8 2H8"), ModificationSites.K);

                    if ((Form1.NEUCODE_SIXPLEX_ARG && countResidues('R', peptide.Sequence) < 1) || Form1.NEUCODE_SIXPLEX_LEU && countResidues('L', peptide.Sequence) < 1)
                    {
                        theoMasses[0, 0] = check1.Mass.Monoisotopic;
                        theoMasses[1, 0] = check2.Mass.Monoisotopic;
                    }
                    else
                    {
                        check3.SetModification(NamedChemicalFormula.GetModification("Lys +8 13C6 15N2"), ModificationSites.K);
                        check4.SetModification(NamedChemicalFormula.GetModification("Lys +8 2H8"), ModificationSites.K);
                        check5.SetModification(NamedChemicalFormula.GetModification("Lys +8 13C6 15N2"), ModificationSites.K);
                        check6.SetModification(NamedChemicalFormula.GetModification("Lys +8 2H8"), ModificationSites.K);

                        if (Form1.NEUCODE_SIXPLEX_ARG)
                        {
                            check3.SetModification(NamedChemicalFormula.GetModification("Arg +6 13C6"), ModificationSites.R);
                            check4.SetModification(NamedChemicalFormula.GetModification("Arg +6 13C6"), ModificationSites.R);
                            check5.SetModification(NamedChemicalFormula.GetModification("Arg +10 13C6 15N4"), ModificationSites.R);
                            check6.SetModification(NamedChemicalFormula.GetModification("Arg +10 13C6 15N4"), ModificationSites.R);
                        }
                        else
                        {
                            check3.SetModification(NamedChemicalFormula.GetModification("Leu +6 13C6 15N1"), ModificationSites.L);
                            check4.SetModification(NamedChemicalFormula.GetModification("Leu +6 13C6 15N1"), ModificationSites.L);
                            check5.SetModification(NamedChemicalFormula.GetModification("Leu +10 2H10"), ModificationSites.L);
                            check6.SetModification(NamedChemicalFormula.GetModification("Leu +10 2H10"), ModificationSites.L);
                        }

                        theoMasses[0, 0] = check1.Mass.Monoisotopic;
                        theoMasses[1, 0] = check2.Mass.Monoisotopic;
                        theoMasses[2, 0] = check3.Mass.Monoisotopic;
                        theoMasses[3, 0] = check4.Mass.Monoisotopic;
                        theoMasses[4, 0] = check5.Mass.Monoisotopic;
                        theoMasses[5, 0] = check6.Mass.Monoisotopic;
                    }
                }

                else if (Form1.SILAC_DUPLEX_LYSCN || Form1.SILAC_DUPLEX_LYSH || Form1.SILAC_DUPLEX_LEUCN || Form1.SILAC_DUPLEX_LEUH)
                {
                    check1 = new Peptide(peptide);
                    check2 = new Peptide(peptide);

                    if (Form1.SILAC_DUPLEX_LYSCN)
                    {
                        check2.SetModification(NamedChemicalFormula.GetModification("Lys +8 13C6 15N2"), ModificationSites.K);
                    }
                    else if (Form1.SILAC_DUPLEX_LYSH)
                    {
                        check2.SetModification(NamedChemicalFormula.GetModification("Lys +8 2H8"), ModificationSites.K);
                    }
                    else if (Form1.SILAC_DUPLEX_LEUCN)
                    {
                        check2.SetModification(NamedChemicalFormula.GetModification("Leu +7 13C6 15N1"), ModificationSites.L);
                    }
                    else
                    {
                        check2.SetModification(NamedChemicalFormula.GetModification("Leu +7 2H7"), ModificationSites.L);
                    }

                    theoMasses[0, 0] = check1.Mass.Monoisotopic;
                    theoMasses[1, 0] = check2.Mass.Monoisotopic;
                }
            }
        }

        /* Counts the number of specific residues in the peptide
         */
        public int countResidues(char residue, string sequence)
        {
            string res = residue.ToString();
            if (res != res.ToUpper()) res = res.ToUpper();
            int count = 0;
            for (int i = 0; i < sequence.Length; i++)
            {
                if (sequence[i].Equals(residue))
                {
                    count++;
                }
            }
            return count;
        }

        public int findLowerResolutionCharge(ILabeledPeak peak, MSDataScan scan)
        {
            int lowerResolutionCharge = 0;
            List<MZPeak> peaks = null;
            SortedList<double, ILabeledPeak> closestPeaks;
            MassTolerance tolerance = new MassTolerance(MassToleranceType.PPM, 10.0);

            scan.MassSpectrum.TryGetPeaks(peak.X - tolerance.Value, peak.X + tolerance.Value, out peaks);

            if (peaks.Count > 0)
            {
                closestPeaks = new SortedList<double, ILabeledPeak>();
                foreach (MZPeak peakToConvert in peaks)
                {
                    ILabeledPeak convertedPeak = (ILabeledPeak)peakToConvert;
                    try
                    {
                        closestPeaks.Add(Math.Abs(convertedPeak.X - peak.X), convertedPeak);
                    }
                    catch (ArgumentException)
                    {
                    }
                }
                lowerResolutionCharge = closestPeaks.ElementAt(0).Value.Charge;
            }

            return lowerResolutionCharge;
        }

        /* Finds and returns the largest peak within the specified m/z range
         */
        public ILabeledPeak largestPeak (double theoMass, MSDataScan currentScan, MassTolerance range, MSDataFile rawFile)
        {
            PeptideSpectralMatch best = bestPSMs[rawFile];
            List<MZPeak> allPeaks;
            List<ILabeledPeak> allConvertedPeaks;
            List<ILabeledPeak> sortedIntensityPeaks;
            ILabeledPeak largest = null;
            int charge = best.Charge;
            double theoMZ = Mass.MzFromMass(theoMass, charge);
            MassRange mZMassRange = new MassRange(theoMZ, range);
            double minSN = 10.0;
            if (Form1.FIRSTSEARCHDONE)
            {
                minSN = Form1.MINIMUMSN;
            }

            if (theoMass == 0)
            {
                return null;
            }
            else
            {
                // Find the largest peak in the range
                if (currentScan.MassSpectrum.TryGetPeaks(mZMassRange, out allPeaks))
                {
                    allConvertedPeaks = new List<ILabeledPeak>();
                    foreach (MZPeak peak in allPeaks)
                    {
                        allConvertedPeaks.Add((ILabeledPeak)peak);
                    }
                    sortedIntensityPeaks = allConvertedPeaks.OrderByDescending(peakSort => peakSort.GetSignalToNoise()).ToList();
                    largest = sortedIntensityPeaks[0];
                }
                // Verify that the largest peak meets the signal-to-noise and charge state criteria
                if (largest != null)
                {
                    int peakCharge = 0;
                    // Ignore any peak with a known charge state that is incorrect
                    if (largest.Charge != charge)
                    {
                        MSDataScan previous = rawFile[currentScan.SpectrumNumber - 1];
                        peakCharge = findLowerResolutionCharge(largest, previous);
                    }
                    
                    if (peakCharge > 0 && peakCharge != charge) 
                    {
                        largest = null;
                    }
                    // Ignore any peak with a signal-to-noise below the specified threshold
                    else if (largest.GetSignalToNoise() < minSN)
                    {
                        largest = null;
                    }
                    else
                    {
                        return largest;
                    }
                }
            }
            return null;
        }

        /* Searches for sets of isotopologues and, when found, calculates the ppm error associated with each isotopologue
         */
        public void patternPPM (double[] theoMasses, MSDataScan currentScan, MassTolerance range, MSDataFile rawFile, List<PrecursorPPM> ppms)
        {
            PeptideSpectralMatch best = bestPSMs[rawFile];
            int charge = best.Charge;
            double theoMass = 0;

            // Find the average theoretical mass for the set of isotopologues
            for (int i = 0; i < theoMasses.Length; i++)
            {
                theoMass += theoMasses[i];
            }

            theoMass = theoMass / ((double)theoMasses.Length);

            MassTolerance tolerance = new MassTolerance(MassToleranceType.PPM, 50.0);

            MassRange minRange = new MassRange(theoMasses[0], tolerance);
            MassRange maxRange = new MassRange(theoMasses[theoMasses.Length - 1], tolerance);

            double theoMZ = Mass.MzFromMass(theoMass, charge);
            double minMZ = Mass.MzFromMass(minRange.Minimum, charge);
            double maxMZ = Mass.MzFromMass(maxRange.Maximum, charge);
            MassRange mZMassRange = new MassRange(minMZ, maxMZ);

            if (theoMass == 0)
            {
                return;
            }
            else
            {
                List<MZPeak> peaks = null;
                List<ILabeledPeak> sortedPeaks = null;
                List<ILabeledPeak> topIntensityPeaks = null;
                List<ILabeledPeak> mzOrderedPeaks = null;
                if (currentScan.MassSpectrum.TryGetPeaks(mZMassRange, out peaks))
                {
                    sortedPeaks = new List<ILabeledPeak>();
                    foreach (IPeak peak in peaks)
                    {
                        ILabeledPeak newPeak = (ILabeledPeak)peak;
                        if (newPeak.GetSignalToNoise() < 10.0)
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
                    topIntensityPeaks = sortedPeaks.OrderByDescending(peak => peak.GetSignalToNoise()).Take(numIsotopologues).ToList();

                    // Only consider complete sets of isotopologues
                    if (topIntensityPeaks.Count == numIsotopologues)
                    {
                        mzOrderedPeaks = topIntensityPeaks.OrderBy(peak => peak.X).ToList();
                        PrecursorPPM ppm;

                        for (int i = 0; i < theoMasses.Length; i++)
                        {
                            ppm = new PrecursorPPM(charge, this.sequence, best.EValue, MassTolerance.GetTolerance(mzOrderedPeaks[i].X, Mass.MzFromMass(theoMasses[i], charge), MassToleranceType.PPM));
                            ppms.Add(ppm);
                        }
                    }
                }
            }
        }

        /* Calculates the MS1 scan range to use for quantification
         * Includes all PSMs, extending above and below according to the entered retention time window
         */
        public void calculateScanRange(MSDataFile rawFile, double rTWindow)
        {
            int best = bestPSMs[rawFile].ScanNumber;
            List<PeptideSpectralMatch> sortedPSMs = PSMs[rawFile].OrderBy(PSM => PSM.ScanNumber).ToList();
            int first = sortedPSMs[0].ScanNumber;
            int last = sortedPSMs[sortedPSMs.Count - 1].ScanNumber;
            double firstRunTime = rawFile.GetRetentionTime(rawFile.FirstSpectrumNumber); //Earliest time in run
            double lastRunTime = rawFile.GetRetentionTime(rawFile.LastSpectrumNumber); //Latest time in run
            double firstScanTime = rawFile.GetRetentionTime(first);
            double lastScanTime = rawFile.GetRetentionTime(last);
            double bestScanTime = rawFile.GetRetentionTime(best);

            fullScanList = new List<MSDataScan>();

            // For peptides eluting across more than 1000 scans, create retention time window around best PSM
            if (last - first < 1000)
            {
                rTRange = new Range<double>(firstScanTime - rTWindow, lastScanTime + rTWindow);
            }
            else
            {
                rTRange = new Range<double>(bestScanTime - rTWindow, bestScanTime + rTWindow);
            }

            // Check to see if retention time window extends past the length of the run
            if (rTRange.Minimum < firstRunTime)
            {
                rTRange.Minimum = firstRunTime; //Minimum time surpassed
            }

            if (rTRange.Maximum > lastRunTime)
            {
                rTRange.Maximum = lastRunTime; //Maximum time surpassed
            }

            // For NeuCode, use high resolution MS1 scans (> 200,000)
            if (Form1.MULTIINJECT)
            {
                foreach (PeptideSpectralMatch psm in sortedPSMs)
                {
                    int correctScan = -1;
                    bool stop = false;
                    int scanCount = psm.ScanNumber + 1;
                    while (scanCount < rawFile.LastSpectrumNumber && !stop)
                    {
                        if (rawFile.GetMsScan(scanCount).MsnOrder == 2 && ((MsnDataScan)rawFile.GetMsScan(scanCount)).DissociationType == DissociationType.HCD)
                        {
                            correctScan = scanCount;
                            fullScanList.Add(rawFile[correctScan]);
                            stop = true;
                        }
                        scanCount++;
                    }
                }
            }
            else if (Form1.NEUCODE || Form1.SILAC_DUPLEX_LEUCN || Form1.SILAC_DUPLEX_LEUH)
            {
                foreach (MSDataScan scan in rawFile.Where(scan => scan.MsnOrder == 1 && scan.Resolution > 200000 && scan.RetentionTime <= rTRange.Maximum && scan.RetentionTime >= rTRange.Minimum))
                {
                    fullScanList.Add(scan);
                }
            }
            //Otherwise, include every MS1 scan for quantification
            else
            {
                foreach (MSDataScan scan in rawFile.Where(scan => scan.MsnOrder == 1 && scan.RetentionTime <= rTRange.Maximum && scan.RetentionTime >= rTRange.Minimum))
                {
                    fullScanList.Add(scan);
                }
            }
        }

        /* Finds light and heavy peaks in the current scan, using a specified PPM window centered around the peptide's adjusted masses
         * Adds a pair to the peptide's list of all pairs for that raw file if enough (i.e., equal to or greater than the number of channels) peaks are found
         */
        public void findPeaks(MSDataScan current, MSDataFile rawFile)
        {
            int charge = bestPSMs[rawFile].Charge;
            ILabeledPeak[,] peaksFound = new ILabeledPeak[numChannels, numIsotopes];
            Pair pair;
            double injectionTime = rawFile.GetInjectionTime(current.SpectrumNumber);

            // Count the number of non-null peaks returned
            int peaksCount = 0;
            int peaksPerIsotope = 0;

            //if (Form1.NEUCODE_ARG_PROLINECONVERSION)
            //{
            //    double prolineConversion = 5.016775;
            //    ILabeledPeak prolinePeak;
            //    for (int j = 0; j < numIsotopes; j++)
            //    {
            //        peaksFound[0, j] = largestPeak(adjustedTheoMasses[0, j], current, tolerance, rawFile);
            //        if (peaksFound[0, j] != null)
            //        {
            //            peaksCount++;
            //            // Look for proline conversion
            //            if (prolineCount > 0)
            //            {
            //                for (int p = 1; p <= prolineCount; p++)
            //                {
            //                    prolinePeak = largestPeak(adjustedTheoMasses[0, j] + (p * prolineConversion), current, tolerance, rawFile);
            //                    if (prolinePeak != null)
            //                    {
            //                        if (peaksFound[1, j] == null)
            //                        {
            //                            peaksFound[1, j] = prolinePeak;
            //                        }
            //                        else
            //                        {
            //                            ThermoLabeledPeak duplicatePeak = (ThermoLabeledPeak)peaksFound[1, j];
            //                            ThermoLabeledPeak addIntensities = new ThermoLabeledPeak(duplicatePeak.MZ, duplicatePeak.Intensity + prolinePeak.Y, duplicatePeak.Charge, duplicatePeak.Noise);
            //                        }
            //                    }
            //                }
            //            }
            //            // Look for "fake" proline conversion in peptides that do not contain proline
            //            else
            //            {
            //                int fakeProlines = 2;
            //                for (int p = 1; p <= fakeProlines; p++)
            //                {
            //                    prolinePeak = largestPeak(adjustedTheoMasses[0, j] + (p * prolineConversion), current, tolerance, rawFile);
            //                    if (prolinePeak != null)
            //                    {
            //                        if (peaksFound[1, j] == null)
            //                        {
            //                            peaksFound[1, j] = prolinePeak;
            //                        }
            //                        else
            //                        {
            //                            ThermoLabeledPeak duplicatePeak = (ThermoLabeledPeak)peaksFound[1, j];
            //                            ThermoLabeledPeak addIntensities = new ThermoLabeledPeak(duplicatePeak.MZ, duplicatePeak.Intensity + prolinePeak.Y, duplicatePeak.Charge, duplicatePeak.Noise);
            //                        }
            //                    }
            //                }
            //            }
            //        }
            //    }

            //    // Create a new pair associated with the MS1 scan and peptide
            //    pair = new Pair(this, rawFile, current.SpectrumNumber);
            //    pair.peaks = peaksFound;
            //    if (peaksCount >= peaksNeeded)
            //    {
            //        allPairs[rawFile].Add(pair);
            //    }
            //}
            if (Form1.NEUCODE)
            {
                List<ILabeledPeak> nonNullPeaks;
                int missingChannelCount;
                int duplicatePeakCount;
                int channelIndex;
                ILabeledPeak peak;
                int clusterCount = 0; // # of clusters meeting the peak count criteria 
                double patternRange = 0.005;

                for (int c = 0; c < numClusters; c++) // Start cluster loop
                {
                    for (int j = 0; j < numIsotopes; j++) // Start isotope loop
                    {
                        missingChannelCount = 0;
                        duplicatePeakCount = 0;
                        nonNullPeaks = new List<ILabeledPeak>();
                        channelIndex = c * numIsotopologues;

                        for (int i = channelIndex; i < channelIndex + numIsotopologues; i++) // Start isotopologue loop
                        {
                            peak = largestPeak(adjustedTheoMasses[i, j], current, tolerance, rawFile);
                            if (peak != null)
                            {
                                if (nonNullPeaks.Contains(peak))
                                {
                                    duplicatePeakCount++;
                                }
                                else
                                {
                                    nonNullPeaks.Add(peak);
                                    peaksFound[i, j] = peak;
                                }
                            }
                            else
                            {
                                missingChannelCount++;
                            }
                        }

                        // Search for patterns if missing or duplicate peaks detected
                        if (duplicatePeakCount > 0 || missingChannelCount > 0)
                        {
                            double min = Mass.MzFromMass(adjustedTheoMasses[channelIndex, j], charge) - patternRange;
                            double max = Mass.MzFromMass(adjustedTheoMasses[channelIndex + (numIsotopologues - 1), j], charge) + patternRange;
                            List<MZPeak> peaksReTry = null;
                            List<ILabeledPeak> convertedPeaks = null;
                            List<ILabeledPeak> topSignalToNoisePeaks = null;
                            List<ILabeledPeak> orderedByMzPeaks = null;
                            if (current.MassSpectrum.TryGetPeaks(min, max, out peaksReTry))
                            {
                                convertedPeaks = new List<ILabeledPeak>();
                                foreach (IPeak peakReTry in peaksReTry)
                                {
                                    ILabeledPeak newPeak = (ILabeledPeak)peakReTry;
                                    if (newPeak.GetSignalToNoise() < SIGNALTONOISE)
                                    {
                                        newPeak = null;
                                    }
                                    else if (newPeak.Charge != 0 && newPeak.Charge != charge)
                                    {
                                        newPeak = null;
                                    }
                                    else
                                    {
                                        convertedPeaks.Add(newPeak);
                                    }
                                }
                                topSignalToNoisePeaks = convertedPeaks.OrderByDescending(peakSort => peakSort.GetSignalToNoise()).Take(numIsotopologues).ToList();

                                // Peaks found for a complete set of isotopologues                             
                                if (topSignalToNoisePeaks.Count == numIsotopologues)
                                {
                                    orderedByMzPeaks = topSignalToNoisePeaks.OrderBy(peakSort => peakSort.X).ToList();
                                    for (int i = 0; i < numIsotopologues; i++)
                                    {
                                        peaksFound[channelIndex + i, j] = orderedByMzPeaks[i];
                                    }
                                }
                                else if (Form1.NOISEBANDCAP)
                                {
                                    // Peaks found for an incomplete set of isotopologues by pattern detection
                                    if (topSignalToNoisePeaks.Count >= peaksNeeded)
                                    {
                                        ILabeledPeak[] mappedPeaks = mapPeaks(topSignalToNoisePeaks, channelIndex, channelIndex + (numIsotopologues - 1), j, charge);
                                        for (int m = channelIndex; m < channelIndex + numIsotopologues; m++)
                                        {
                                            peaksFound[m, j] = mappedPeaks[m - channelIndex];
                                        }
                                    }
                                    // Peaks found for an incomplete set of isotopologues by peak detection
                                    else if (nonNullPeaks.Count >= peaksNeeded)
                                    {
                                        ILabeledPeak[] mappedPeaks = mapPeaks(nonNullPeaks, channelIndex, channelIndex + (numIsotopologues - 1), j, charge);
                                        for (int m = channelIndex; m < channelIndex + numIsotopologues; m++)
                                        {
                                            peaksFound[m, j] = mappedPeaks[m - channelIndex];
                                        }
                                    }
                                }
                                else
                                {
                                    for (int i = 0; i < numIsotopologues; i++)
                                    {
                                        peaksFound[channelIndex + i, j] = null;
                                    }
                                }
                            }
                        } // End pattern detection

                        // Count detected peaks
                        peaksCount = 0;
                        for (int i = channelIndex; i < channelIndex + numIsotopologues; i++)
                        {
                            if (peaksFound[i, j] != null)
                            {
                                peaksCount++;
                            }
                        }

                        // If enough peaks are detected, increase the cluster count
                        if (peaksCount >= peaksNeeded)
                        {
                            clusterCount++;
                        }
                        // If not, set all peaks to null
                        else
                        {
                            for (int i = channelIndex; i < channelIndex + numIsotopologues; i++)
                            {
                                peaksFound[i, j] = null;
                            }
                        }
                    } // End isotope loop
                } // End cluster loop

                if (clusterCount > 0)
                {
                    pair = new Pair(this, rawFile, current.SpectrumNumber);
                    pair.peaks = peaksFound;
                    pair.injectionTime = injectionTime;
                    allPairs[rawFile].Add(pair);
                }   
            }
            else
            {
                for (int j = 0; j < numIsotopes; j++)
                {
                    peaksPerIsotope = 0;
                    for (int i = 0; i < numChannels; i++)
                    {
                        peaksFound[i, j] = largestPeak(adjustedTheoMasses[i, j], current, tolerance, rawFile);
                        if (peaksFound[i, j] != null)
                        {
                            peaksPerIsotope++;
                        }
                    }
                    if (peaksPerIsotope < numChannels && !Form1.NOISEBANDCAP)
                    {
                        for (int k = 0; k < numChannels; k++)
                        {
                            peaksFound[k, j] = null;
                        }
                    }
                }

                for (int i = 0; i < numChannels; i++)
                {
                    for (int j = 0; j < numIsotopes; j++)
                    {
                        if (peaksFound[i, j] != null)
                        {
                            peaksCount++;
                        }
                    }
                }

                if (peaksCount >= peaksNeeded)
                {
                    pair = new Pair(this, rawFile, current.SpectrumNumber);
                    pair.peaks = peaksFound;
                    pair.injectionTime = injectionTime;
                    allPairs[rawFile].Add(pair);
                }
            }
        }

        /* When given a list of peaks corresponding to an incomplete set of isotopologues, finds the combination with the smallest associated ppm error
         */
        public ILabeledPeak[] mapPeaks(List<ILabeledPeak> peaks, int channelStart, int channelEnd, int isotope, int charge)
        {
            int numChannels = channelEnd - channelStart + 1;
            ILabeledPeak[] mappedPeaks = new ILabeledPeak[numChannels];
            List<ILabeledPeak> orderedByMzPeaks = null;

            orderedByMzPeaks = peaks.OrderBy(peak => peak.X).ToList();

            int possibleCombinations = 0;
            if (peaks.Count == 1)
            {
                if (numChannels == 2) possibleCombinations = 2;
                else if (numChannels == 3) possibleCombinations = 3;
            }
            else if (peaks.Count == 2)
            {
                if (numChannels == 3) possibleCombinations = 3;
                else if (numChannels == 4) possibleCombinations = 6;
            }
            else if (peaks.Count == 3)
            {
                if (numChannels == 4) possibleCombinations = 4;
                else if (numChannels == 6) possibleCombinations = 20;
            }
            else if (peaks.Count == 4) possibleCombinations = 14;
            else if (peaks.Count == 5) possibleCombinations = 6;
            KeyValuePair<double, ILabeledPeak>[,] errorCheck = new KeyValuePair<double, ILabeledPeak>[possibleCombinations, numChannels];
            
            // Then calculate all the errors
            ILabeledPeak currentPeak;
            double neutralMass;
            double error1;
            double error2;
            double error3;
            double error4;
            double error5;
            double error6;
            double overallError;
            double minError;
            int minErrorIndex;
            for (int i = 0; i < peaks.Count; i++)
            {
                currentPeak = peaks[i];

                if (numIsotopologues == 2)
                {
                    neutralMass = Mass.MassFromMz(currentPeak.X, charge);
                    error1 = MassTolerance.GetTolerance(neutralMass, adjustedTheoMasses[channelStart, isotope], MassToleranceType.PPM);
                    error2 = MassTolerance.GetTolerance(neutralMass, adjustedTheoMasses[channelEnd, isotope], MassToleranceType.PPM);
                    if (Math.Abs(error1) < Math.Abs(error2))
                    {
                        mappedPeaks[0] = currentPeak;
                    }
                    else
                    {
                        mappedPeaks[1] = currentPeak;
                    }
                    return mappedPeaks;
                }
                else if (numIsotopologues == 3)
                {
                    neutralMass = Mass.MassFromMz(currentPeak.X, charge);
                    error1 = MassTolerance.GetTolerance(neutralMass, adjustedTheoMasses[channelStart, isotope], MassToleranceType.PPM);
                    error2 = MassTolerance.GetTolerance(neutralMass, adjustedTheoMasses[channelStart + 1, isotope], MassToleranceType.PPM);
                    error3 = MassTolerance.GetTolerance(neutralMass, adjustedTheoMasses[channelStart + 2, isotope], MassToleranceType.PPM);

                    if (peaks.Count == 1)
                    {
                        if (error1 < error2 && error1 < error3) mappedPeaks[0] = currentPeak;
                        else if (error2 < error1 && error2 < error3) mappedPeaks[1] = currentPeak;
                        else if (error3 < error1 && error3 < error2) mappedPeaks[2] = currentPeak;
                    }
                    else if (peaks.Count == 2)
                    {
                        if (i == 0)
                        {
                            errorCheck[0, 0] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error1), currentPeak);
                            errorCheck[1, 0] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error1), currentPeak);
                            errorCheck[2, 1] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error2), currentPeak);
                        }
                        else if (i == 1)
                        {
                            errorCheck[0, 1] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error2), currentPeak);
                            errorCheck[1, 2] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error3), currentPeak);
                            errorCheck[2, 2] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error3), currentPeak);
                        }
                    }
                }
                else if (numIsotopologues == 4)
                {
                    neutralMass = Mass.MassFromMz(currentPeak.X, charge);
                    error1 = MassTolerance.GetTolerance(neutralMass, adjustedTheoMasses[channelStart, isotope], MassToleranceType.PPM);
                    error2 = MassTolerance.GetTolerance(neutralMass, adjustedTheoMasses[channelStart + 1, isotope], MassToleranceType.PPM);
                    error3 = MassTolerance.GetTolerance(neutralMass, adjustedTheoMasses[channelStart + 2, isotope], MassToleranceType.PPM);
                    error4 = MassTolerance.GetTolerance(neutralMass, adjustedTheoMasses[channelStart + 3, isotope], MassToleranceType.PPM);

                    if (peaks.Count == 2)
                    {
                        if (i == 0)
                        {
                            errorCheck[0, 0] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error1), currentPeak);
                            errorCheck[1, 0] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error1), currentPeak);
                            errorCheck[2, 0] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error1), currentPeak);
                            errorCheck[3, 1] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error2), currentPeak);
                            errorCheck[4, 1] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error2), currentPeak);
                            errorCheck[5, 2] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error3), currentPeak);
                        }
                        else
                        {
                            errorCheck[0, 1] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error2), currentPeak);
                            errorCheck[1, 2] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error3), currentPeak);
                            errorCheck[2, 3] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error4), currentPeak);
                            errorCheck[3, 2] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error3), currentPeak);
                            errorCheck[4, 3] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error4), currentPeak);
                            errorCheck[5, 3] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error4), currentPeak);
                        }
                    }

                    if (peaks.Count == 3)
                    {
                        if (i == 0)
                        {
                            errorCheck[0, 0] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error1), currentPeak);
                            errorCheck[1, 0] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error1), currentPeak);
                            errorCheck[2, 0] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error1), currentPeak);
                            errorCheck[3, 1] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error2), currentPeak);
                        }
                        else if (i == 1)
                        {
                            errorCheck[0, 1] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error2), currentPeak);
                            errorCheck[1, 1] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error2), currentPeak);
                            errorCheck[2, 2] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error3), currentPeak);
                            errorCheck[3, 2] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error3), currentPeak);
                        }
                        else
                        {
                            errorCheck[0, 2] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error3), currentPeak);
                            errorCheck[1, 3] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error4), currentPeak);
                            errorCheck[2, 3] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error4), currentPeak);
                            errorCheck[3, 3] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error4), currentPeak);
                        }
                    }
                }
                else if (numIsotopologues == 6)
                {
                    neutralMass = Mass.MassFromMz(currentPeak.X, charge);
                    error1 = MassTolerance.GetTolerance(neutralMass, adjustedTheoMasses[channelStart, isotope], MassToleranceType.PPM);
                    error2 = MassTolerance.GetTolerance(neutralMass, adjustedTheoMasses[channelStart + 1, isotope], MassToleranceType.PPM);
                    error3 = MassTolerance.GetTolerance(neutralMass, adjustedTheoMasses[channelStart + 2, isotope], MassToleranceType.PPM);
                    error4 = MassTolerance.GetTolerance(neutralMass, adjustedTheoMasses[channelStart + 3, isotope], MassToleranceType.PPM);
                    error5 = MassTolerance.GetTolerance(neutralMass, adjustedTheoMasses[channelStart + 4, isotope], MassToleranceType.PPM);
                    error6 = MassTolerance.GetTolerance(neutralMass, adjustedTheoMasses[channelStart + 5, isotope], MassToleranceType.PPM);

                    if (peaks.Count == 3)
                    {
                        if (i == 0)
                        {
                            errorCheck[0, 0] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error1), currentPeak);
                            errorCheck[1, 0] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error1), currentPeak);
                            errorCheck[2, 0] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error1), currentPeak);
                            errorCheck[3, 0] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error1), currentPeak);
                            errorCheck[4, 0] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error1), currentPeak);
                            errorCheck[5, 0] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error1), currentPeak);
                            errorCheck[6, 0] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error1), currentPeak);
                            errorCheck[7, 0] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error1), currentPeak);
                            errorCheck[8, 0] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error1), currentPeak);
                            errorCheck[9, 0] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error1), currentPeak);
                            errorCheck[10, 1] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error2), currentPeak);
                            errorCheck[11, 1] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error2), currentPeak);
                            errorCheck[12, 1] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error2), currentPeak);
                            errorCheck[13, 1] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error2), currentPeak);
                            errorCheck[14, 1] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error2), currentPeak);
                            errorCheck[15, 1] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error2), currentPeak);
                            errorCheck[16, 2] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error3), currentPeak);
                            errorCheck[17, 2] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error3), currentPeak);
                            errorCheck[18, 2] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error3), currentPeak);
                            errorCheck[19, 3] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error4), currentPeak);
                            
                        }
                        else if (i == 1)
                        {
                            errorCheck[0, 1] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error2), currentPeak);
                            errorCheck[1, 1] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error2), currentPeak);
                            errorCheck[2, 1] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error2), currentPeak);
                            errorCheck[3, 1] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error2), currentPeak);
                            errorCheck[4, 2] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error3), currentPeak);
                            errorCheck[5, 2] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error3), currentPeak);
                            errorCheck[6, 2] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error3), currentPeak);
                            errorCheck[7, 3] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error4), currentPeak);
                            errorCheck[8, 3] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error4), currentPeak);
                            errorCheck[9, 4] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error5), currentPeak);
                            errorCheck[10, 2] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error3), currentPeak);
                            errorCheck[11, 2] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error3), currentPeak);
                            errorCheck[12, 2] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error3), currentPeak);
                            errorCheck[13, 3] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error4), currentPeak);
                            errorCheck[14, 3] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error4), currentPeak);
                            errorCheck[15, 4] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error5), currentPeak);
                            errorCheck[16, 3] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error4), currentPeak);
                            errorCheck[17, 3] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error4), currentPeak);
                            errorCheck[18, 4] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error5), currentPeak);
                            errorCheck[19, 4] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error5), currentPeak);
                        }
                        else if (i == 2)
                        {
                            errorCheck[0, 2] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error3), currentPeak);
                            errorCheck[1, 3] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error4), currentPeak);
                            errorCheck[2, 4] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error5), currentPeak);
                            errorCheck[3, 5] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error6), currentPeak);
                            errorCheck[4, 3] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error4), currentPeak);
                            errorCheck[5, 4] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error5), currentPeak);
                            errorCheck[6, 5] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error6), currentPeak);
                            errorCheck[7, 4] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error5), currentPeak);
                            errorCheck[8, 5] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error6), currentPeak);
                            errorCheck[9, 5] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error6), currentPeak);
                            errorCheck[10, 3] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error4), currentPeak);
                            errorCheck[11, 4] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error5), currentPeak);
                            errorCheck[12, 5] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error6), currentPeak);
                            errorCheck[13, 4] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error5), currentPeak);
                            errorCheck[14, 5] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error6), currentPeak);
                            errorCheck[15, 5] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error6), currentPeak);
                            errorCheck[16, 4] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error5), currentPeak);
                            errorCheck[17, 5] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error6), currentPeak);
                            errorCheck[18, 5] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error6), currentPeak);
                            errorCheck[19, 5] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error6), currentPeak);
                        }
                    }
                    else if (peaks.Count == 4)
                    {
                        if (i == 0)
                        {
                            errorCheck[0, 0] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error1), currentPeak);
                            errorCheck[1, 0] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error1), currentPeak);
                            errorCheck[2, 0] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error1), currentPeak);
                            errorCheck[3, 0] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error1), currentPeak);
                            errorCheck[4, 0] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error1), currentPeak);
                            errorCheck[5, 0] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error2), currentPeak);
                            errorCheck[6, 0] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error1), currentPeak);
                            errorCheck[7, 0] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error1), currentPeak);
                            errorCheck[8, 0] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error1), currentPeak);
                            errorCheck[9, 1] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error2), currentPeak);
                            errorCheck[10, 1] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error2), currentPeak);
                            errorCheck[11, 1] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error2), currentPeak);
                            errorCheck[12, 1] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error2), currentPeak);
                            errorCheck[13, 2] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error3), currentPeak);
                        }
                        else if (i == 1)
                        {
                            errorCheck[0, 1] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error2), currentPeak);
                            errorCheck[1, 1] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error2), currentPeak);
                            errorCheck[2, 1] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error2), currentPeak);
                            errorCheck[3, 1] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error2), currentPeak);
                            errorCheck[4, 1] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error2), currentPeak);
                            errorCheck[5, 1] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error2), currentPeak);
                            errorCheck[6, 2] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error3), currentPeak);
                            errorCheck[7, 2] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error3), currentPeak);
                            errorCheck[8, 3] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error4), currentPeak);
                            errorCheck[9, 2] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error3), currentPeak);
                            errorCheck[10, 2] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error3), currentPeak);
                            errorCheck[11, 2] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error3), currentPeak);
                            errorCheck[12, 3] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error4), currentPeak);
                            errorCheck[13, 3] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error4), currentPeak);
                        }
                        else if (i == 2)
                        {
                            errorCheck[0, 2] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error3), currentPeak);
                            errorCheck[1, 2] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error3), currentPeak);
                            errorCheck[2, 2] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error3), currentPeak);
                            errorCheck[3, 3] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error4), currentPeak);
                            errorCheck[4, 3] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error4), currentPeak);
                            errorCheck[5, 4] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error5), currentPeak);
                            errorCheck[6, 3] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error4), currentPeak);
                            errorCheck[7, 3] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error4), currentPeak);
                            errorCheck[8, 4] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error5), currentPeak);
                            errorCheck[9, 3] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error4), currentPeak);
                            errorCheck[10, 3] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error4), currentPeak);
                            errorCheck[11, 4] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error5), currentPeak);
                            errorCheck[12, 4] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error5), currentPeak);
                            errorCheck[13, 4] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error5), currentPeak);
                        }
                        else if (i == 3)
                        {
                            errorCheck[0, 3] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error4), currentPeak);
                            errorCheck[1, 4] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error5), currentPeak);
                            errorCheck[2, 5] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error6), currentPeak);
                            errorCheck[3, 4] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error5), currentPeak);
                            errorCheck[4, 5] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error6), currentPeak);
                            errorCheck[5, 5] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error6), currentPeak);
                            errorCheck[6, 4] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error5), currentPeak);
                            errorCheck[7, 5] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error6), currentPeak);
                            errorCheck[8, 5] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error6), currentPeak);
                            errorCheck[9, 4] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error5), currentPeak);
                            errorCheck[10, 5] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error6), currentPeak);
                            errorCheck[11, 5] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error6), currentPeak);
                            errorCheck[12, 5] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error6), currentPeak);
                            errorCheck[13, 5] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error6), currentPeak);
                        }
                    }
                    else if (peaks.Count == 5)
                    {
                        if (i == 0)
                        {
                            errorCheck[0, 0] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error1), currentPeak);
                            errorCheck[1, 0] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error1), currentPeak);
                            errorCheck[2, 0] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error1), currentPeak);
                            errorCheck[3, 0] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error1), currentPeak);
                            errorCheck[4, 0] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error1), currentPeak);
                            errorCheck[5, 1] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error2), currentPeak);
                        }
                        else if (i == 1)
                        {
                            errorCheck[0, 1] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error2), currentPeak);
                            errorCheck[1, 1] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error2), currentPeak);
                            errorCheck[2, 1] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error2), currentPeak);
                            errorCheck[3, 1] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error2), currentPeak);
                            errorCheck[4, 2] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error3), currentPeak);
                            errorCheck[5, 2] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error3), currentPeak);
                        }
                        else if (i == 2)
                        {
                            errorCheck[0, 2] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error3), currentPeak);
                            errorCheck[1, 2] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error3), currentPeak);
                            errorCheck[2, 2] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error3), currentPeak);
                            errorCheck[3, 3] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error4), currentPeak);
                            errorCheck[4, 3] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error4), currentPeak);
                            errorCheck[5, 3] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error4), currentPeak);
                        }
                        else if (i == 3)
                        {
                            errorCheck[0, 3] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error4), currentPeak);
                            errorCheck[1, 3] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error4), currentPeak);
                            errorCheck[2, 4] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error5), currentPeak);
                            errorCheck[3, 4] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error5), currentPeak);
                            errorCheck[4, 4] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error5), currentPeak);
                            errorCheck[5, 4] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error5), currentPeak);
                        }
                        else if (i == 4)
                        {
                            errorCheck[0, 4] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error5), currentPeak);
                            errorCheck[1, 5] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error6), currentPeak);
                            errorCheck[2, 5] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error6), currentPeak);
                            errorCheck[3, 5] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error6), currentPeak);
                            errorCheck[4, 5] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error6), currentPeak);
                            errorCheck[5, 5] = new KeyValuePair<double, ILabeledPeak>(Math.Abs(error6), currentPeak);
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

            // Set mappedPeaks to peak combination with lowest overall error
            for (int n = 0; n < numChannels; n++)
            {
                mappedPeaks[n] = errorCheck[minErrorIndex, n].Value;
            }

            return mappedPeaks;
        }

        /* Checks a pair for coalescence by considering the missing channel frequencies of more and less intense pairs
         */
        public void checkPairCoalescence(MSDataFile rawFile)
        {
            PeptideSpectralMatch best = PSMs[rawFile][0];
            int scanNumber = best.ScanNumber;
            double injectionTime;
            List<Pair> pairs;
            Pair pair;
            ILabeledPeak peak1;
            ILabeledPeak peak2;
            ILabeledPeak peak3;
            ILabeledPeak peak4;
            IsotopePair isotopePair;
            List<IsotopePair> lowerIntensity;
            List<IsotopePair> higherIntensity;
            double lowerFrequency;
            double higherFrequency;

            allPairs.TryGetValue(rawFile, out pairs);
            for (int p = 0; p < pairs.Count(); p++)
            {
                pair = pairs[p];
                injectionTime = rawFile.GetInjectionTime(pair.scanNumber);
                for (int c = 0; c < numChannels; c += numIsotopologues)
                {
                    for (int j = 0; j < numIsotopes; j++)
                    {
                        if (numIsotopologues == 2)
                        {
                            peak1 = pair.peaks[c, j];
                            peak2 = pair.peaks[c + 1, j];

                            // Do not alter complete or null pairs or pairs with a maximum intensity below the set threshold
                            if (peak1 != null && peak2 != null)
                            {

                            }
                            else if (peak1 == null && peak2 == null)
                            {
                                
                            }
                            else if (maximumIntensity <= System.Math.Log10(Form1.MAXIMUMDNL))
                            {
                                
                            }
                            // Proceed with the analysis for all other pairs
                            else
                            {
                                ILabeledPeak singlePeak;
                                ILabeledPeak otherPeak1;
                                ILabeledPeak otherPeak2;
                                double singlePeakIntensity;
                                Pair otherPair;

                                // Find the single peak
                                if (peak1 == null)
                                {
                                    singlePeak = peak2;
                                }
                                else
                                {
                                    singlePeak = peak1;
                                }
                                singlePeakIntensity = singlePeak.GetDenormalizedIntensity(injectionTime);

                                // Sort the peptide's other pairs (complete and incomplete) based on their intensities relative to the single peak
                                lowerIntensity = new List<IsotopePair>();
                                higherIntensity = new List<IsotopePair>();
                                for (int a = 0; a < pairs.Count(); a++)
                                {
                                    otherPair = pairs[a];
                                    for (int h = 0; h < numChannels; h += numIsotopologues)
                                    {
                                        for (int s = 0; s < numIsotopes; s++)
                                        {
                                            otherPeak1 = otherPair.peaks[h, s];
                                            otherPeak2 = otherPair.peaks[h + 1, s];

                                            //Both channels present
                                            if (otherPeak1 != null && otherPeak2 != null)
                                            {
                                                isotopePair = new IsotopePair(otherPair, s, otherPeak1.X, otherPeak1.GetDenormalizedIntensity(injectionTime), otherPeak2.X, otherPeak2.GetDenormalizedIntensity(injectionTime));
                                                if (isotopePair.intensity <= singlePeakIntensity)
                                                {
                                                    lowerIntensity.Add(isotopePair);
                                                }
                                                else
                                                {
                                                    higherIntensity.Add(isotopePair);
                                                }
                                            }
                                            else if (otherPeak1 == null && otherPeak2 == null)
                                            {

                                            }
                                            //One channel present
                                            else
                                            {
                                                if (otherPeak1 != null)
                                                {
                                                    isotopePair = new IsotopePair(otherPair, s, otherPeak1.X, otherPeak1.GetDenormalizedIntensity(injectionTime), 0, 0);
                                                    if (isotopePair.intensity <= singlePeakIntensity)
                                                    {
                                                        lowerIntensity.Add(isotopePair);
                                                    }
                                                    else
                                                    {
                                                        higherIntensity.Add(isotopePair);
                                                    }
                                                }
                                                if (otherPeak2 != null)
                                                {
                                                    isotopePair = new IsotopePair(otherPair, s, otherPeak2.X, otherPeak2.GetDenormalizedIntensity(injectionTime), 0, 0);
                                                    if (isotopePair.intensity <= singlePeakIntensity)
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

                                // Calculate the missing channel frequencies for pairs both lower and higher in intensity than the single peak
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

                                // A peak is deemed coalesced if the frequency of a missing channel is at least 1.5-fold greater for the more intense pairs than the less intense (each category must have at least 2 contributing pairs)
                                if (lowerIntensity.Count() > 1 && higherIntensity.Count() > 1 && higherFrequency / lowerFrequency > 1.5)
                                {
                                    coalescenceDetected = true;
                                    if (coalescedPeakIntensities == null)
                                    {
                                        coalescedPeakIntensities = new List<double>();
                                    }
                                    else
                                    {
                                        // Keep track of peak intensities that are deemed to be coalesced
                                        coalescedPeakIntensities.Add(singlePeakIntensity);
                                    }

                                    // Set single peak to null
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

                            // Do nothing for complete pairs, null pairs, or pairs below the intensity threshold
                            if (peak1 != null && peak2 != null && peak3 != null && peak4 != null)
                            {
                                
                            }
                            else if (peak1 == null && peak2 == null && peak3 == null && peak4 == null)
                            {
                                
                            }
                            else if (maximumIntensity <= System.Math.Log10(Form1.MAXIMUMDNL))
                            {
                                
                            }
                            // Proceed with all other pairs
                            else 
                            {
                                ILabeledPeak otherPeak1;
                                ILabeledPeak otherPeak2;
                                ILabeledPeak otherPeak3;
                                ILabeledPeak otherPeak4;
                                double presentPeaksTotalIntensity = 0;
                                int presentPeaksCount = 0;
                                double presentPeaksIntensity = 0;

                                // Calculate the maximum intensity of all present peaks
                                if (peak1 != null)
                                {
                                    if (peak1.GetDenormalizedIntensity(injectionTime) > presentPeaksIntensity)
                                    {
                                        presentPeaksTotalIntensity += peak1.GetDenormalizedIntensity(injectionTime);
                                        presentPeaksCount++;
                                    }
                                }
                                if (peak2 != null)
                                {
                                    if (peak2.GetDenormalizedIntensity(injectionTime) > presentPeaksIntensity)
                                    {
                                        presentPeaksTotalIntensity += peak2.GetDenormalizedIntensity(injectionTime);
                                        presentPeaksCount++;
                                    }
                                }
                                if (peak3 != null)
                                {
                                    if (peak3.GetDenormalizedIntensity(injectionTime) > presentPeaksIntensity)
                                    {
                                        presentPeaksTotalIntensity = peak3.GetDenormalizedIntensity(injectionTime);
                                        presentPeaksCount++;
                                    }
                                }
                                if (peak4 != null)
                                {
                                    if (peak4.GetDenormalizedIntensity(injectionTime) > presentPeaksIntensity)
                                    {
                                        presentPeaksTotalIntensity = peak4.GetDenormalizedIntensity(injectionTime);
                                        presentPeaksCount++;
                                    }
                                }

                                presentPeaksIntensity = presentPeaksTotalIntensity / ((double)presentPeaksCount);

                                // Sort the peptide's other pairs (either both or one channel present) based on their intensities relative to the single peak
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
                                            otherPeak1 = otherPair.peaks[h, s];
                                            otherPeak2 = otherPair.peaks[h + 1, s];
                                            otherPeak3 = otherPair.peaks[h + 2, s];
                                            otherPeak4 = otherPair.peaks[h + 3, s];

                                            // All channels present
                                            if (otherPeak1 != null && otherPeak2 != null && otherPeak3 != null && otherPeak4 != null)
                                            {
                                                isotopePair = new IsotopePair(otherPair, s, otherPeak1.X, otherPeak1.GetDenormalizedIntensity(injectionTime), otherPeak2.X, otherPeak2.GetDenormalizedIntensity(injectionTime), otherPeak3.X, otherPeak3.GetDenormalizedIntensity(injectionTime), otherPeak4.X, otherPeak4.GetDenormalizedIntensity(injectionTime));
                                                if (isotopePair.intensity <= presentPeaksIntensity)
                                                {
                                                    lowerIntensity.Add(isotopePair);
                                                }
                                                else
                                                {
                                                    higherIntensity.Add(isotopePair);
                                                }
                                            }
                                            // No channels present
                                            else if (otherPeak1 == null && otherPeak2 == null && otherPeak3 == null && otherPeak4 == null)
                                            {
                                                // Do nothing
                                            }
                                            // Incomplete pairs -- set null peaks to blank peaks
                                            else
                                            {
                                                MZPeak blank = new MZPeak(0,0);

                                                if (otherPeak1 == null)
                                                {
                                                    otherPeak1 = (ILabeledPeak)blank;
                                                }
                                                if (otherPeak2 == null)
                                                {
                                                    otherPeak2 = (ILabeledPeak)blank;
                                                }
                                                if (otherPeak3 == null)
                                                {
                                                    otherPeak3 = (ILabeledPeak)blank;
                                                }
                                                if (otherPeak4 == null)
                                                {
                                                    otherPeak4 = (ILabeledPeak)blank;
                                                }

                                                isotopePair = new IsotopePair(otherPair, s, otherPeak1.X, otherPeak1.GetDenormalizedIntensity(injectionTime), otherPeak2.X, otherPeak2.GetDenormalizedIntensity(injectionTime), otherPeak3.X, otherPeak3.GetDenormalizedIntensity(injectionTime), otherPeak4.X, otherPeak4.GetDenormalizedIntensity(injectionTime));
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

                                // Calculate the missing channel frequencies for pairs both lower and higher in intensity than the single peak
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

                                // A pair is deemed coalesced if the frequency of a missing channel is at least 1.5-fold greater for the more intense pairs than the less intense (each category must have at least 2 contributing pairs)
                                if (lowerIntensity.Count() > 1 && higherIntensity.Count() > 1 && higherFrequency / lowerFrequency > 1.5)
                                {
                                    coalescenceDetected = true;
                                    if (coalescedPeakIntensities == null)
                                    {
                                        coalescedPeakIntensities = new List<double>();
                                    }
                                    else
                                    {
                                        coalescedPeakIntensities.Add(presentPeaksIntensity);
                                    }

                                    // Set single peaks to null for a coalesced pair
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

        /* Checks the pair spacing to make sure it falls within the spacing tolerance
         */
        public void checkPairSpacing(MSDataFile rawFile, List<Spacing> spacings)
        {
            int charge = bestPSMs[rawFile].Charge;
            List<Pair> pairs;
            Pair pair;
            ILabeledPeak peak1;
            ILabeledPeak peak2;
            ILabeledPeak peak3;
            ILabeledPeak peak4;
            ILabeledPeak peak5;
            ILabeledPeak peak6;
            double spacing1;
            double spacing2;
            double spacing3;
            double spacing4;
            double spacing5;
            double peak1NeutralMass;
            double peak2NeutralMass;
            double peak3NeutralMass;
            double peak4NeutralMass;
            double peak5NeutralMass;
            double peak6NeutralMass;

            allPairs.TryGetValue(rawFile, out pairs);
            if (numIsotopologues == 1 || numIsotopologues == 2)
            {
                for (int i = 0; i < pairs.Count; i++)
                {
                    pair = pairs[i];
                    for (int j = 0; j < numChannels; j += 2)
                    {
                        for (int k = 0; k < numIsotopes; k++)
                        {
                            peak1 = pair.peaks[j, k];
                            peak2 = pair.peaks[j + 1, k];

                            // Only look at complete pairs
                            if (peak1 != null && peak2 != null)
                            {
                                peak1NeutralMass = Mass.MassFromMz(peak1.X, charge);
                                peak2NeutralMass = Mass.MassFromMz(peak2.X, charge);

                                spacing1 = peak2NeutralMass - peak1NeutralMass;
                                //theoSpacing = (double)lysineCount * lysineVersion;
                                //Spacing space = new Spacing(theoSpacing, spacing, charge, light.MZ, heavy.MZ);
                                //spacings.Add(space);
                                //Remove pairs whose spacing is more than 5 mDa away from the calculated spacing
                                if (spacing1 < spacingMassRange[0].Minimum || spacing1 > spacingMassRange[0].Maximum)
                                {
                                    //Console.WriteLine("bad spacing");
                                    pair.peaks[j, k] = null;
                                    pair.peaks[j + 1, k] = null;
                                }
                            }
                            //else if ((peak1 != null || peak2 != null) && Form1.NOISEBANDCAP)
                            //{
                            //    // Do nothing
                            //}
                            else
                            {
                                pair.peaks[j, k] = null;
                                pair.peaks[j + 1, k] = null;
                            }
                        }
                    }
                }
            }
            else if (numIsotopologues == 3)
            {
                for (int i = 0; i < pairs.Count; i++)
                {
                    pair = pairs[i];
                    for (int j = 0; j < numChannels; j += 3)
                    {
                        for (int k = 0; k < numIsotopes; k++)
                        {
                            peak1 = pair.peaks[j, k];
                            peak2 = pair.peaks[j + 1, k];
                            peak3 = pair.peaks[j + 2, k];

                            // Only look at complete pairs
                            if (peak1 != null && peak2 != null && peak3 != null)
                            {
                                peak1NeutralMass = Mass.MassFromMz(peak1.X, charge);
                                peak2NeutralMass = Mass.MassFromMz(peak2.X, charge);
                                peak3NeutralMass = Mass.MassFromMz(peak3.X, charge);

                                spacing1 = peak2NeutralMass - peak1NeutralMass;
                                spacing2 = peak3NeutralMass - peak2NeutralMass;
                                //theoSpacing = (double)lysineCount * lysineVersion;
                                //Spacing space = new Spacing(theoSpacing, spacing, charge, light.MZ, heavy.MZ);
                                //spacings.Add(space);
                                //Remove pairs whose spacing is more than 5 mDa away from the calculated spacing
                                if (spacing1 < spacingMassRange[0].Minimum || spacing1 > spacingMassRange[0].Maximum)
                                {
                                    //Console.WriteLine("bad spacing");
                                    pair.peaks[j, k] = null;
                                    if (!Form1.NOISEBANDCAP)
                                    {
                                        pair.peaks[j + 1, k] = null;
                                        pair.peaks[j + 2, k] = null;
                                    }
                                }
                                if (spacing2 < spacingMassRange[1].Minimum || spacing2 > spacingMassRange[1].Maximum)
                                {
                                    pair.peaks[j + 2, k] = null;
                                    if (pair.peaks[j, k] == null) pair.peaks[j + 1, k] = null;

                                    if (!Form1.NOISEBANDCAP)
                                    {
                                        pair.peaks[j, k] = null;
                                        pair.peaks[j + 1, k] = null;
                                    }
                                }
                            }
                            //else if ((peak1 != null || peak2 != null) && Form1.NOISEBANDCAP)
                            //{
                            //    // Do nothing
                            //}
                            else
                            {
                                pair.peaks[j, k] = null;
                                pair.peaks[j + 1, k] = null;
                                pair.peaks[j + 2, k] = null;
                            }
                        }
                    }
                }
            }
            else if (numIsotopologues == 4)
            {
                for (int i = 0; i < pairs.Count; i++)
                {
                    pair = pairs[i];
                    for (int j = 0; j < numChannels; j += 4)
                    {
                        for (int k = 0; k < numIsotopes; k++)
                        {
                            peak1 = pair.peaks[j, k];
                            peak2 = pair.peaks[j + 1, k];
                            peak3 = pair.peaks[j + 2, k];
                            peak4 = pair.peaks[j + 3, k];

                            // Only look at complete pairs
                            if (peak1 != null && peak2 != null && peak3 != null)
                            {
                                peak1NeutralMass = Mass.MassFromMz(peak1.X, charge);
                                peak2NeutralMass = Mass.MassFromMz(peak2.X, charge);
                                peak3NeutralMass = Mass.MassFromMz(peak3.X, charge);
                                peak4NeutralMass = Mass.MassFromMz(peak4.X, charge);

                                spacing1 = peak2NeutralMass - peak1NeutralMass;
                                spacing2 = peak3NeutralMass - peak2NeutralMass;
                                spacing3 = peak4NeutralMass - peak3NeutralMass;
                                //theoSpacing = (double)lysineCount * lysineVersion;
                                //Spacing space = new Spacing(theoSpacing, spacing, charge, light.MZ, heavy.MZ);
                                //spacings.Add(space);
                                //Remove pairs whose spacing is more than 5 mDa away from the calculated spacing
                                if (spacing1 < spacingMassRange[0].Minimum || spacing1 > spacingMassRange[0].Maximum)
                                {
                                    //Console.WriteLine("bad spacing");
                                    pair.peaks[j, k] = null;
                                    if (!Form1.NOISEBANDCAP)
                                    {
                                        pair.peaks[j + 1, k] = null;
                                        pair.peaks[j + 2, k] = null;
                                        pair.peaks[j + 3, k] = null;
                                    }
                                }
                                if (spacing3 < spacingMassRange[2].Minimum || spacing3 > spacingMassRange[2].Maximum)
                                {
                                    pair.peaks[j + 3, k] = null;
                                    if (!Form1.NOISEBANDCAP)
                                    {
                                        pair.peaks[j, k] = null;
                                        pair.peaks[j + 1, k] = null;
                                        pair.peaks[j + 2, k] = null;
                                    }
                                }
                                if (spacing2 < spacingMassRange[1].Minimum || spacing2 > spacingMassRange[1].Maximum)
                                {
                                    if (pair.peaks[j, k] == null) pair.peaks[j + 1, k] = null;
                                    if (pair.peaks[j + 3, k] == null) pair.peaks[j + 2, k] = null;

                                    if (!Form1.NOISEBANDCAP)
                                    {
                                        pair.peaks[j, k] = null;
                                        pair.peaks[j + 1, k] = null;
                                        pair.peaks[j + 2, k] = null;
                                        pair.peaks[j + 3, k] = null;
                                    }
                                }
                            }
                            //else if ((peak1 != null || peak2 != null) && Form1.NOISEBANDCAP)
                            //{
                            //    // Do nothing
                            //}
                            else
                            {
                                pair.peaks[j, k] = null;
                                pair.peaks[j + 1, k] = null;
                                pair.peaks[j + 2, k] = null;
                                pair.peaks[j + 3, k] = null;
                            }
                        }
                    }
                }
            }
            else if (numIsotopologues == 6)
            {
                for (int i = 0; i < pairs.Count; i++)
                {
                    pair = pairs[i];
                    for (int j = 0; j < numChannels; j += 6)
                    {
                        for (int k = 0; k < numIsotopes; k++)
                        {
                            peak1 = pair.peaks[j, k];
                            peak2 = pair.peaks[j + 1, k];
                            peak3 = pair.peaks[j + 2, k];
                            peak4 = pair.peaks[j + 3, k];
                            peak5 = pair.peaks[j + 4, k];
                            peak6 = pair.peaks[j + 5, k];

                            // Only look at complete pairs
                            if (peak1 != null && peak2 != null && peak3 != null && peak4 != null && peak5 != null && peak6 != null)
                            {
                                peak1NeutralMass = Mass.MassFromMz(peak1.X, charge);
                                peak2NeutralMass = Mass.MassFromMz(peak2.X, charge);
                                peak3NeutralMass = Mass.MassFromMz(peak3.X, charge);
                                peak4NeutralMass = Mass.MassFromMz(peak4.X, charge);
                                peak5NeutralMass = Mass.MassFromMz(peak5.X, charge);
                                peak6NeutralMass = Mass.MassFromMz(peak6.X, charge);

                                spacing1 = peak2NeutralMass - peak1NeutralMass;
                                spacing2 = peak3NeutralMass - peak2NeutralMass;
                                spacing3 = peak4NeutralMass - peak3NeutralMass;
                                spacing4 = peak5NeutralMass - peak4NeutralMass;
                                spacing5 = peak6NeutralMass - peak5NeutralMass;
                                //theoSpacing = (double)lysineCount * lysineVersion;
                                //Spacing space = new Spacing(theoSpacing, spacing, charge, light.MZ, heavy.MZ);
                                //spacings.Add(space);
                                //Remove pairs whose spacing is more than 5 mDa away from the calculated spacing
                                if (spacing1 < spacingMassRange[0].Minimum || spacing1 > spacingMassRange[0].Maximum)
                                {
                                    //Console.WriteLine("bad spacing");
                                    pair.peaks[j, k] = null;
                                    if (!Form1.NOISEBANDCAP)
                                    {
                                        pair.peaks[j + 1, k] = null;
                                        pair.peaks[j + 2, k] = null;
                                        pair.peaks[j + 3, k] = null;
                                        pair.peaks[j + 4, k] = null;
                                        pair.peaks[j + 5, k] = null;
                                    }
                                }
                                if (spacing5 < spacingMassRange[4].Minimum || spacing5 > spacingMassRange[4].Maximum)
                                {
                                    pair.peaks[j + 5, k] = null;
                                    if (!Form1.NOISEBANDCAP)
                                    {
                                        pair.peaks[j, k] = null;
                                        pair.peaks[j + 1, k] = null;
                                        pair.peaks[j + 2, k] = null;
                                        pair.peaks[j + 3, k] = null;
                                        pair.peaks[j + 4, k] = null;
                                    }
                                }
                                if (spacing2 < spacingMassRange[1].Minimum || spacing2 > spacingMassRange[1].Maximum)
                                {
                                    if (pair.peaks[j, k] == null) pair.peaks[j + 1, k] = null;
                                    if (!Form1.NOISEBANDCAP)
                                    {
                                        pair.peaks[j, k] = null;
                                        pair.peaks[j + 1, k] = null;
                                        pair.peaks[j + 2, k] = null;
                                        pair.peaks[j + 3, k] = null;
                                        pair.peaks[j + 4, k] = null;
                                        pair.peaks[j + 5, k] = null;
                                    }
                                }
                                if (spacing4 < spacingMassRange[3].Minimum || spacing4 > spacingMassRange[3].Maximum)
                                {
                                    if (pair.peaks[j + 5, k] == null) pair.peaks[j + 4, k] = null;
                                    if (!Form1.NOISEBANDCAP)
                                    {
                                        pair.peaks[j, k] = null;
                                        pair.peaks[j + 1, k] = null;
                                        pair.peaks[j + 2, k] = null;
                                        pair.peaks[j + 3, k] = null;
                                        pair.peaks[j + 4, k] = null;
                                        pair.peaks[j + 5, k] = null;
                                    }
                                }
                                if (spacing3 < spacingMassRange[2].Minimum || spacing3 > spacingMassRange[2].Maximum)
                                {
                                    if (pair.peaks[j + 1, k] == null) pair.peaks[j + 2, k] = null;
                                    if (pair.peaks[j + 4, k] == null) pair.peaks[j + 3, k] = null;
                                    if (!Form1.NOISEBANDCAP)
                                    {
                                        pair.peaks[j, k] = null;
                                        pair.peaks[j + 1, k] = null;
                                        pair.peaks[j + 2, k] = null;
                                        pair.peaks[j + 3, k] = null;
                                        pair.peaks[j + 4, k] = null;
                                        pair.peaks[j + 5, k] = null;
                                    }
                                }
                            }
                            //else if ((peak1 != null || peak2 != null) && Form1.NOISEBANDCAP)
                            //{
                            //    // Do nothing
                            //}
                            else
                            {
                                pair.peaks[j, k] = null;
                                pair.peaks[j + 1, k] = null;
                                pair.peaks[j + 2, k] = null;
                                pair.peaks[j + 3, k] = null;
                                pair.peaks[j + 4, k] = null;
                                pair.peaks[j + 5, k] = null;
                            }
                        }
                    }
                }
            }
            if (!Form1.NOISEBANDCAP)
            {
                List<Pair> completePairList;
                completePairs.TryGetValue(rawFile, out completePairList);

                foreach (Pair currentPair in pairs)
                {
                    if (currentPair.totalPeakCount > 0)
                    {
                        completePairList.Add(currentPair);
                    }
                }
                pairs.Clear();

            }

            

            //if (numIsotopologues == 4)
            //{
            //    for (int p = 0; p < pairs.Count(); p++)
            //    {
            //        pair = pairs[p];
            //        int peakCount;

            //        for (int c = 0; c < numClusters; c++)
            //        {
            //            for (int j = 0; j < numIsotopes; j ++)
            //            {
            //                peakCount = 0;
            //                int channelIndex = c * numIsotopologues;
            //                for (int i = c; i < c + numIsotopologues; i += numIsotopologues)
            //                { 
            //                    peak1 = pair.peaks[channelIndex, j];
            //                    peak2 = pair.peaks[channelIndex + 1, j];
            //                    peak3 = pair.peaks[channelIndex + 2, j];
            //                    peak4 = pair.peaks[channelIndex + 3, j];
            //                    if (peak1 != null)
            //                    {
            //                        peakCount++;
            //                    }
            //                    if (peak2 != null)
            //                    {
            //                        peakCount++;
            //                    }
            //                    if (peak3 != null)
            //                    {
            //                        peakCount++;
            //                    }
            //                    if (peak4 != null)
            //                    {
            //                        peakCount++;
            //                    }

            //                    // First consider complete pairs
            //                    if (peakCount == 4)
            //                    {
            //                        peak1NeutralMass = Mass.MassFromMz(peak1.X, charge);
            //                        peak2NeutralMass = Mass.MassFromMz(peak2.X, charge);
            //                        peak3NeutralMass = Mass.MassFromMz(peak3.X, charge);
            //                        peak4NeutralMass = Mass.MassFromMz(peak4.X, charge);

            //                        spacing1 = peak4NeutralMass - peak3NeutralMass;
            //                        spacing2 = peak3NeutralMass - peak2NeutralMass;
            //                        spacing3 = peak2NeutralMass - peak1NeutralMass;

            //                        if (Form1.NOISEBANDCAP)
            //                        {
            //                            if (spacing1 < spacingMassRange.Minimum || spacing1 > spacingMassRange.Maximum)
            //                            {
            //                                pair.peaks[channelIndex + 3, j] = null;
            //                                peakCount--;
            //                            }
            //                            if (spacing2 < spacingMassRange.Minimum || spacing2 > spacingMassRange.Maximum)
            //                            {
            //                                pair.peaks[channelIndex + 2, j] = null;
            //                                peakCount--;
            //                            }
            //                            if (spacing3 < spacingMassRange.Minimum || spacing3 > spacingMassRange.Maximum)
            //                            {
            //                                pair.peaks[channelIndex, j] = null;
            //                                peakCount--;
            //                                pair.peaks[channelIndex + 1, j] = null;
            //                                peakCount--;
            //                            }

            //                            // Verify that enough peaks are still present
            //                            if (peakCount < peaksNeeded)
            //                            {
            //                                pair.peaks[channelIndex, j] = null;
            //                                pair.peaks[channelIndex + 1, j] = null;
            //                                pair.peaks[channelIndex + 2, j] = null;
            //                                pair.peaks[channelIndex + 3, j] = null;
            //                            }
            //                        }
            //                        else
            //                        {
            //                            if (spacing1 < spacingMassRange.Minimum || spacing1 > spacingMassRange.Maximum || spacing1 < spacingMassRange.Minimum || spacing1 > spacingMassRange.Maximum || spacing3 < spacingMassRange.Minimum || spacing3 > spacingMassRange.Maximum)
            //                            {
            //                                pair.peaks[channelIndex, j] = null;
            //                                pair.peaks[channelIndex + 1, j] = null;
            //                                pair.peaks[channelIndex + 2, j] = null;
            //                                pair.peaks[channelIndex + 3, j] = null;
            //                            }
            //                        }
            //                    }
            //                    else if (peakCount == 3)
            //                    {
            //                        if (Form1.NOISEBANDCAP)
            //                        {
            //                            if (peak1 != null)
            //                            {
            //                                peak1NeutralMass = Mass.MassFromMz(peak1.X, charge);
            //                                if (peak2 != null && peak3 != null) // peaks 1,2,3
            //                                {
            //                                    peak2NeutralMass = Mass.MassFromMz(peak2.X, charge);
            //                                    peak3NeutralMass = Mass.MassFromMz(peak3.X, charge);
            //                                    spacing1 = peak3NeutralMass - peak2NeutralMass;
            //                                    spacing2 = peak2NeutralMass - peak1NeutralMass;
            //                                    if (spacing1 > spacingMassRange.Maximum || spacing1 < spacingMassRange.Minimum)
            //                                    {
            //                                        pair.peaks[channelIndex + 2, j] = null;
            //                                        peakCount--;
            //                                    }
            //                                    if (spacing2 > spacingMassRange.Maximum || spacing2 < spacingMassRange.Minimum)
            //                                    {
            //                                        pair.peaks[channelIndex + 1, j] = null;
            //                                        peakCount--;
            //                                        pair.peaks[channelIndex, j] = null;
            //                                        peakCount--;
            //                                    }
            //                                }
            //                                else if (peak3 != null && peak4 != null) // peaks 1,3,4
            //                                {
            //                                    peak3NeutralMass = Mass.MassFromMz(peak3.X, charge);
            //                                    peak4NeutralMass = Mass.MassFromMz(peak4.X, charge);
            //                                    spacing1 = peak4NeutralMass - peak3NeutralMass;
            //                                    spacing2 = (peak3NeutralMass - peak1NeutralMass) / 2.0;
            //                                    if (spacing1 > spacingMassRange.Maximum || spacing1 < spacingMassRange.Minimum)
            //                                    {
            //                                        pair.peaks[channelIndex + 3, j] = null;
            //                                        peakCount--;
            //                                    }
            //                                    if (spacing2 > spacingMassRange.Maximum || spacing2 < spacingMassRange.Minimum)
            //                                    {
            //                                        pair.peaks[channelIndex + 2, j] = null;
            //                                        peakCount--;
            //                                        pair.peaks[channelIndex, j] = null;
            //                                        peakCount--;
            //                                    }
            //                                }
            //                                else // peaks 1,2,4
            //                                {
            //                                    peak2NeutralMass = Mass.MassFromMz(peak2.X, charge);
            //                                    peak4NeutralMass = Mass.MassFromMz(peak4.X, charge);
            //                                    spacing1 = (peak4NeutralMass - peak2NeutralMass) / 2.0;
            //                                    spacing2 = peak2NeutralMass - peak1NeutralMass;
            //                                    if (spacing1 > spacingMassRange.Maximum || spacing1 < spacingMassRange.Minimum)
            //                                    {
            //                                        pair.peaks[channelIndex + 3, j] = null;
            //                                        peakCount--;
            //                                    }
            //                                    if (spacing2 > spacingMassRange.Maximum || spacing2 < spacingMassRange.Minimum)
            //                                    {
            //                                        pair.peaks[channelIndex + 1, j] = null;
            //                                        peakCount--;
            //                                        pair.peaks[channelIndex, j] = null;
            //                                        peakCount--;
            //                                    }
            //                                }
            //                            }
            //                            else 
            //                            {
            //                                peak2NeutralMass = Mass.MassFromMz(peak2.X, charge);
            //                                peak3NeutralMass = Mass.MassFromMz(peak3.X, charge);
            //                                peak4NeutralMass = Mass.MassFromMz(peak4.X, charge);
            //                                spacing1 = peak4NeutralMass - peak3NeutralMass;
            //                                spacing2 = peak3NeutralMass - peak2NeutralMass;
            //                                if (spacing1 > spacingMassRange.Maximum || spacing1 < spacingMassRange.Minimum)
            //                                {
            //                                    pair.peaks[channelIndex + 3, j] = null;
            //                                    peakCount--;
            //                                }
            //                                if (spacing2 > spacingMassRange.Maximum || spacing2 < spacingMassRange.Minimum)
            //                                {
            //                                    pair.peaks[channelIndex + 2, j] = null;
            //                                    peakCount--;
            //                                    pair.peaks[channelIndex + 1, j] = null;
            //                                    peakCount--;
            //                                }
            //                            }

            //                            // Verify that enough peaks are still present
            //                            if (peakCount < peaksNeeded)
            //                            {
            //                                pair.peaks[channelIndex, j] = null;
            //                                pair.peaks[channelIndex + 1, j] = null;
            //                                pair.peaks[channelIndex + 2, j] = null;
            //                                pair.peaks[channelIndex + 3, j] = null;
            //                            }
            //                        }
            //                        else
            //                        {
            //                            pair.peaks[channelIndex, j] = null;
            //                            pair.peaks[channelIndex + 1, j] = null;
            //                            pair.peaks[channelIndex + 2, j] = null;
            //                            pair.peaks[channelIndex + 3, j] = null;
            //                        }
            //                    }
            //                    else if (peakCount == 2)
            //                    {
            //                        if (Form1.NOISEBANDCAP)
            //                        {
            //                            if (peak1 != null)
            //                            {
            //                                peak1NeutralMass = Mass.MassFromMz(peak1.X, charge);
            //                                if (peak2 != null)
            //                                {
            //                                    peak2NeutralMass = Mass.MassFromMz(peak2.X, charge);
            //                                    spacing1 = peak2NeutralMass - peak1NeutralMass;
            //                                }
            //                                else if (peak3 != null)
            //                                {
            //                                    peak3NeutralMass = Mass.MassFromMz(peak3.X, charge);
            //                                    spacing1 = (peak3NeutralMass - peak1NeutralMass) / 2.0;
            //                                }
            //                                else
            //                                {
            //                                    peak4NeutralMass = Mass.MassFromMz(peak4.X, charge);
            //                                    spacing1 = (peak4NeutralMass - peak1NeutralMass) / 3.0;
            //                                }
            //                            }
            //                            else if (peak2 != null)
            //                            {
            //                                peak2NeutralMass = Mass.MassFromMz(peak2.X, charge);
            //                                if (peak3 != null)
            //                                {
            //                                    peak3NeutralMass = Mass.MassFromMz(peak3.X, charge);
            //                                    spacing1 = peak3NeutralMass - peak2NeutralMass;
            //                                }
            //                                else if (peak4 != null)
            //                                {
            //                                    peak4NeutralMass = Mass.MassFromMz(peak4.X, charge);
            //                                    spacing1 = (peak4NeutralMass - peak2NeutralMass) / 2.0;
            //                                }
            //                                else
            //                                {
            //                                    peak1NeutralMass = Mass.MassFromMz(peak1.X, charge);
            //                                    spacing1 = peak2NeutralMass - peak1NeutralMass;
            //                                }
            //                            }
            //                            else if (peak3 != null)
            //                            {
            //                                peak3NeutralMass = Mass.MassFromMz(peak3.X, charge);
            //                                if (peak4 != null)
            //                                {
            //                                    peak4NeutralMass = Mass.MassFromMz(peak4.X, charge);
            //                                    spacing1 = peak4NeutralMass - peak3NeutralMass;
            //                                }
            //                                else if (peak1 != null)
            //                                {
            //                                    peak1NeutralMass = Mass.MassFromMz(peak1.X, charge);
            //                                    spacing1 = (peak3NeutralMass - peak1NeutralMass) / 2.0;
            //                                }
            //                                else
            //                                {
            //                                    peak2NeutralMass = Mass.MassFromMz(peak2.X, charge);
            //                                    spacing1 = peak3NeutralMass - peak2NeutralMass;
            //                                }
            //                            }
            //                            else
            //                            {
            //                                peak4NeutralMass = Mass.MassFromMz(peak4.X, charge);
            //                                if (peak1 != null)
            //                                {
            //                                    peak1NeutralMass = Mass.MassFromMz(peak1.X, charge);
            //                                    spacing1 = (peak4NeutralMass - peak1NeutralMass) / 3.0;
            //                                }
            //                                else if (peak2 != null)
            //                                {
            //                                    peak2NeutralMass = Mass.MassFromMz(peak2.X, charge);
            //                                    spacing1 = (peak4NeutralMass - peak2NeutralMass) / 2.0;
            //                                }
            //                                else
            //                                {
            //                                    peak3NeutralMass = Mass.MassFromMz(peak3.X, charge);
            //                                    spacing1 = peak4NeutralMass - peak3NeutralMass;
            //                                }
            //                            }
            //                            if (spacing1 > spacingMassRange.Maximum || spacing1 < spacingMassRange.Minimum)
            //                            {
            //                                pair.peaks[channelIndex, j] = null;
            //                                pair.peaks[channelIndex + 1, j] = null;
            //                                pair.peaks[channelIndex + 2, j] = null;
            //                                pair.peaks[channelIndex + 3, j] = null;
            //                            }
            //                        }
            //                        else
            //                        {
            //                            pair.peaks[channelIndex, j] = null;
            //                            pair.peaks[channelIndex + 1, j] = null;
            //                            pair.peaks[channelIndex + 2, j] = null;
            //                            pair.peaks[channelIndex + 3, j] = null;
            //                        }
            //                    }
            //                    else // For 0 or 1 peaks, no further action needed
            //                    {
            //                        pair.peaks[channelIndex, j] = null;
            //                        pair.peaks[channelIndex + 1, j] = null;
            //                        pair.peaks[channelIndex + 2, j] = null;
            //                        pair.peaks[channelIndex + 3, j] = null;
            //                    }
                    //        } // End isotopologue loop
                    //    } // End isotope loop
                    //} // End cluster loop
                //}
            //}
        }

        /* To missing channels of pairs that are not coalesced, an noise-based intensity is added to permit quantification
         */
        public void applyNoise(MSDataFile rawFile)
        {
            //Apply noise to missing channels
            List<Pair> all;
            List<Pair> completeOnly;
            ILabeledPeak peak;
            double coalescenceIntensity = 100000000.0;
            int channelIndex;

            // If the peptide is found to be coalesced, set the peptide's coalescence intensity threshold at the minimum intensity at which a coalesced pair was detected
            if (coalescenceDetected)
            {
                coalescedPeakIntensities.Sort();
                coalescenceIntensity = coalescedPeakIntensities[0];
            }

            allPairs.TryGetValue(rawFile, out all);
            completePairs.TryGetValue(rawFile, out completeOnly);

            foreach (Pair pair in all)
            {
                for (int j = 0; j < numIsotopes; j++)
                {
                    for (int c = 0; c < numClusters; c++)
                    {
                        channelIndex = c * numIsotopologues;
                        // For complete pairs
                        if (pair.complete[c, j])
                        {
                            // If the pair falls above the peptide's coalescence threshold, set peaks to null
                            if (coalescenceDetected && pair.maxIntensity[c, j] > coalescenceIntensity)
                            {
                                for (int m = channelIndex; m < channelIndex + numIsotopologues; m++)
                                {
                                    pair.peaks[m, j] = null;
                                }
                            }
                            // Otherwise, add pair to complete list
                            else
                            {
                                // First check to see if the peptide is already on the complete list
                                bool found = false;
                                int SN = pair.scanNumber;

                                foreach (Pair noNBC in completeOnly)
                                {
                                    if (noNBC.scanNumber == SN)
                                    {
                                        // If found, update peaks for that pair
                                        for (int m = channelIndex; m < channelIndex + numIsotopologues; m++)
                                        {
                                            noNBC.peaks[m, j] = pair.peaks[m, j];
                                        }
                                        found = true;
                                    }
                                }

                                // If not found, add pair
                                if (!found)
                                {
                                    Pair noNBCPair = new Pair(this, rawFile, pair.scanNumber);
                                    for (int m = channelIndex; m < channelIndex + numIsotopologues; m++)
                                    {
                                        noNBCPair.peaks[m, j] = pair.peaks[m, j];
                                    }
                                    completeOnly.Add(noNBCPair);
                                }
                            }
                        }
                        // For incomplete pairs
                        else
                        {
                            // Check for coalescence
                            if (coalescenceDetected && pair.maxIntensity[c, j] > coalescenceIntensity)
                            {
                                if (pair.maxIntensity[c, j] > coalescenceIntensity)
                                {
                                    for (int m = channelIndex; m < channelIndex + numIsotopologues; m++)
                                    {
                                        pair.peaks[m, j] = null;
                                    }
                                }
                            }
                            // Apply noise to missing channels
                            else
                            {
                                for (int m = channelIndex; m < channelIndex + numIsotopologues; m++)
                                {
                                    if (pair.peaks[m, j] == null)
                                    {
                                        ThermoLabeledPeak noisePeak = new ThermoLabeledPeak(0.0, pair.averageNoise, pair.charge, pair.averageNoise);
                                        pair.peaks[m, j] = noisePeak;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        /* Calculates the systematic error (PPM) associated with the data set
         * Only considers the monoisotopic peaks
         * All isotopologues must be found in order to contribute this calculation
         */
        public void precursorPPMError(MSDataFile rawFile, List<PrecursorPPM> ppms)
        {
            PeptideSpectralMatch best = bestPSMs[rawFile];
            int charge = best.Charge;
            int scanNumber = best.ScanNumber;
            MSDataScan precursorScan;
            ILabeledPeak peak;
            double precursorPPM;
            PrecursorPPM ppm;
            int correctScan = -1;

            // First find the correct precursor scan (i.e., led to the best-scoring PSM)
            if (Form1.MULTIINJECT)
            {
                bool stop = false;
                int scanCount = scanNumber + 1;
                while (scanCount < rawFile.LastSpectrumNumber && !stop)
                {
                    if (rawFile.GetMsScan(scanCount).MsnOrder == 2 && ((MsnDataScan)rawFile.GetMsScan(scanCount)).DissociationType == DissociationType.HCD)
                    {
                        correctScan = scanCount;
                        stop = true;
                    }
                    scanCount++;
                }

                precursorScan = rawFile[correctScan];
            }
            else
            {
                bool stop = false;
                int scanCount = scanNumber - 1;
                while (scanCount < rawFile.LastSpectrumNumber && !stop)
                {
                    if (rawFile.GetMsScan(scanCount).MsnOrder == 1)
                    {
                        correctScan = scanCount;
                        stop = true;
                    }
                    scanCount--;
                }

                precursorScan = rawFile[correctScan];
            }
                
            // Then find all possible precursor monoisotope peaks 

            if (numIsotopologues > 1)
            {
                double[] theo;
                int channelIndex;
                int count;
                MSDataScan lowerResolutionScan = rawFile[correctScan - 1];
                if (Form1.NEUCODE_SIXPLEX_ARG && Form1.CORRECTARGPROLINE && countResidues('R', sequence) > 0 && countResidues('P', sequence) > 0)
                {
                    int correction = correctIsotopeDistributions(lowerResolutionScan, rawFile);
                    theoMasses[4, 0] += correction * (5 * (Constants.CARBON13 - Constants.CARBON));
                    theoMasses[5, 0] += correction * (5 * (Constants.CARBON13 - Constants.CARBON));
                }
                if (Form1.NEUCODE_SIXPLEX_LEU && Form1.CORRECTLEUDLOSS && countResidues('L', sequence) > 0)
                {
                    int correction = correctIsotopeDistributions(lowerResolutionScan, rawFile);
                    theoMasses[4, 0] -= correction * (Constants.DEUTERIUM - Constants.HYDROGEN);
                    theoMasses[5, 0] -= correction * (Constants.DEUTERIUM - Constants.HYDROGEN);
                }
                if (Form1.NEUCODE_DUPLEX_LEU7_18MDA && Form1.CORRECTLEUNLOSS && countResidues('L', sequence) > 0)
                {
                    int correction = correctIsotopeDistributions(lowerResolutionScan, rawFile);
                    theoMasses[0, 0] -= correction * (Constants.NITROGEN15 - Constants.NITROGEN);
                }
                for (int c = 0; c < numClusters; c++)
                {
                    theo = new double[numIsotopologues];
                    channelIndex = c * numIsotopologues;
                    count = 0;
                    for (int i = channelIndex; i < channelIndex + numIsotopologues; i++)
                    {
                        theo[count] = theoMasses[i, 0];
                        count++;
                    }
                    patternPPM(theo, precursorScan, firstSearchMassRange, rawFile, ppms);
                }
            }
            else
            {
                if (Form1.SILAC_DUPLEX_LEUCN && Form1.CORRECTLEUNLOSS)
                {
                    int correction = correctIsotopeDistributions(precursorScan, rawFile);
                    theoMasses[1, 0] -= correction * (Constants.NITROGEN15 - Constants.NITROGEN);
                }
                for (int i = 0; i < numChannels; i++)
                {
                    peak = largestPeak(theoMasses[i, 0], precursorScan, firstSearchMassRange, rawFile);
                    if (peak != null)
                    {
                        precursorPPM = MassTolerance.GetTolerance(Mass.MassFromMz(peak.X, charge), theoMasses[i, 0], MassToleranceType.PPM);
                        ppm = new PrecursorPPM(charge, this.sequence, best.EValue, precursorPPM);
                        ppms.Add(ppm);
                    }
                }
            }             
        }

        public int correctIsotopeDistributions(MSDataScan current, MSDataFile rawFile)
        {
            int numLabels;
            if (Form1.SILAC_DUPLEX_LEUCN || Form1.NEUCODE_DUPLEX_LEU7_18MDA || Form1.NEUCODE_SIXPLEX_LEU)
            {
                numLabels = countResidues('L', sequence) + 1;
            }
            else if (Form1.NEUCODE_SIXPLEX_ARG && countResidues('R', sequence) > 0)
            {
                numLabels = countResidues('P', sequence) + 1;
            }
            else
            {
                numLabels = 0;
            }

            double[,] correctedMasses = new double[numLabels, numIsotopes];// Rows shift theoretical masses to check for isotope conversion, columns are isotopes
            double[,] correctedIntensities = new double[numLabels, numIsotopes]; // Rows shift theoretical masses to check for isotope conversion, columns are isotopes
            int bestMonoMass = 0; // 0 = uncorrected, 1 = corrected for 1 label loss, etc.
            double uncorrectedMonoMass;
            if ((Form1.SILAC_DUPLEX_LEUCN || Form1.NEUCODE_DUPLEX_LEU7_18MDA) && Form1.CORRECTLEUNLOSS && countResidues('L', sequence) > 0)
            {
                // Fill first row
                if (Form1.SILAC_DUPLEX_LEUCN)
                {
                    uncorrectedMonoMass = theoMasses[1, 0];
                }
                else if (Form1.NEUCODE_DUPLEX_LEU7_18MDA)
                {
                    uncorrectedMonoMass = (theoMasses[0, 0] + theoMasses[1, 0]) / 2.0;
                }
                else
                {
                    uncorrectedMonoMass = theoMasses[0, 0];
                }

                // Fill rest of the array
                double mass;
                for (int i = 0; i < numLabels; i++)
                {
                    mass = uncorrectedMonoMass - (i * (Constants.NITROGEN15 - Constants.NITROGEN));
                    for (int j = 0; j < numIsotopes; j++)
                    {
                        correctedMasses[i, j] = mass + (j * (Constants.CARBON13 - Constants.CARBON));
                    }
                }
            }
            else if (Form1.NEUCODE_SIXPLEX_LEU && Form1.CORRECTLEUDLOSS && countResidues('L', sequence) > 0)
            {
                // Fill first row
                uncorrectedMonoMass = (theoMasses[4, 0] + theoMasses[5, 0]) / 2.0;                

                // Fill rest of the array
                double mass;
                for (int i = 0; i < numLabels; i++)
                {
                    mass = uncorrectedMonoMass - (i * (Constants.DEUTERIUM - Constants.HYDROGEN));
                    for (int j = 0; j < numIsotopes; j++)
                    {
                        correctedMasses[i, j] = mass + (j * (Constants.CARBON13 - Constants.CARBON));
                    }
                }
            }
            else if (Form1.NEUCODE_SIXPLEX_ARG && Form1.CORRECTARGPROLINE && countResidues('R', sequence) > 0 && countResidues('P', sequence) > 0)
            {
                double prolineConversionMass = 5 * (Constants.CARBON13 - Constants.CARBON);

                // Fill first row
                uncorrectedMonoMass = (theoMasses[4, 0] + theoMasses[5, 0]) / 2.0;

                // Fill rest of the array
                double mass;
                for (int i = 0; i < numLabels; i++)
                {
                    mass = uncorrectedMonoMass + (i * prolineConversionMass);
                    for (int j = 0; j < numIsotopes; j++)
                    {
                        correctedMasses[i, j] = mass + (j * (Constants.CARBON13 - Constants.CARBON));
                    }
                }
            }

            // Find largest peaks
            ILabeledPeak peak;
            for (int i = 0; i < numLabels; i++)
            {
                for (int j = 0; j < numIsotopes; j++)
                {
                    peak = largestPeak(correctedMasses[i, j], current, firstSearchMassRange, rawFile);
                    if (peak != null)
                    {
                        correctedIntensities[i, j] = peak.Y;
                    }
                }
            }

            // Find best peak combination
            double bestSummedIntensity = 0;
            int bestIndex = 0;
            double sum;
            for (int i = 0; i < numLabels; i++)
            {
                sum = 0;
                for (int j = 0; j < numIsotopes; j++)
                {
                    sum += correctedIntensities[i, j];
                }
                if (sum > bestSummedIntensity)
                {
                    bestSummedIntensity = sum;
                    bestIndex = i;
                }
            }
            bestMonoMass = bestIndex;
            conversionFactor = bestMonoMass;

            return bestMonoMass;
        }

        /* Calculates the proportion of pairs in a given list that have missing channels detected
         * A complete pair has all isotopologues in a cluster present
         * An incomplete pair is missing isotopologues, but must have at least half of the isotopologues detected
         * Example: if 2 possible isotopologues, 1 isotopologue detected is an incomplete pair; if 4 possible isotopologues, 2 or 3 isotopologues detected constitutes a complete pair
         */
        public double calculateMissingChannelFrequency(List<Pair> pairs)
        {
            double frequency = -1;
            int pairTotalCount = 0;
            int incompletePairCount = 0;

            if (pairs != null && pairs.Count > 0)
            {
                foreach (Pair pair in pairs)
                {
                    for (int c = 0; c < numClusters; c++)
                    {
                        for (int j = 0; j < numIsotopes; j++)
                        {
                            // Complete pairs
                            if (pair.peakCount[c, j] == numIsotopologues)
                            {
                                pairTotalCount++;
                            }
                            // Pairs with fewer than half the total number of isotopologues
                            else if (pair.peakCount[c, j] < (numIsotopologues / 2))
                            {
                                // Do nothing
                            }
                            // Incomplete pairs
                            else
                            {
                                pairTotalCount++;
                                incompletePairCount++;
                            }
                        }
                    }
                }
                frequency = ((double)incompletePairCount) / ((double)pairTotalCount);
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

        /* Calculates whether the given pair exceeds the intensity threshold for quantitative filtering
         * Returns true if pair should be included for quantification, false if it should be excluded
         */
        public bool quantFilter(Pair pair, int isotope, int cluster, bool complete)
        {
            bool addIntensities = true;
            double[,] max;
            int channelIndex;

            // Should not have any null peaks in the pair at this point
            if (pair.peakCount[cluster, isotope] != numIsotopologues)
            {
                return false;
            }
            
            // For complete pairs
            if (complete)
            {
                max = maxCompleteIntensity;
                channelIndex = cluster * numIsotopologues;
                for (int i = channelIndex; i < channelIndex + numIsotopologues; i++)
                {
                    double intensityThreshold = INTENSITYCUTOFF * max[i, 0];
                    double denormalizedIntensity = pair.peaks[i, isotope].GetDenormalizedIntensity(pair.injectionTime);
                    if (denormalizedIntensity < intensityThreshold)
                    {
                        return false;
                    }
                }
            }
            // For incomplete pairs
            else
            {
                // Use maximum of complete pairs if more than 2 clusters are complete
                max = maxIntensity;
                channelIndex = cluster * numIsotopologues;
                for (int i = channelIndex; i < channelIndex + numIsotopologues; i++)
                {
                    // For complete pairs
                    if (pair.complete[cluster, isotope])
                    {
                        if (pair.peaks[i, isotope].GetDenormalizedIntensity(pair.injectionTime) < INTENSITYCUTOFF * max[i, 0])
                        {
                            return false;
                        }
                    }
                    // For incomplete pairs, only consider non noise-band capped channels for quantitative filtering
                    if (!pair.complete[cluster, isotope])
                    {
                        if (pair.peaks[i, isotope].X > 0 && pair.peaks[i, isotope].GetDenormalizedIntensity(pair.injectionTime) < INTENSITYCUTOFF * max[i, 0])
                        {
                            return false;
                        }
                    }
                }
            }
            return addIntensities;
        }

        /* Assembles all the peptide's quantitative information
         */
        public void quantify()
        {
            int[,] final = new int[numClusters, 2];
            quantifiedNoiseIncluded = new bool[numClusters];
            finalQuantified = new int[numClusters];
            
            if (completePairs != null && completePairs.Count > 0)
            {
                foreach (List<Pair> pairs in completePairs.Values)
                {
                    foreach (Pair pair in pairs)
                    {
                        for (int j = 0; j < numIsotopes; j++)
                        {
                            for (int c = 0; c < numClusters; c++)
                            {
                                if (Form1.QUANTFILTER)
                                {
                                    // Eliminate low-level pairs by quantitative filtering
                                    if (quantFilter(pair, j, c, true))
                                    {
                                        int channelIndex = c * numIsotopologues;
                                        for (int i = channelIndex; i < channelIndex + numIsotopologues; i++)
                                        {
                                            double denormalizedIntensity = pair.peaks[i, j].GetDenormalizedIntensity(pair.injectionTime);
                                            completeTotalIntensity[i, j] += denormalizedIntensity;
                                            completeTotalIntensity[i, numIsotopes] += denormalizedIntensity;
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
                                            double denormalizedIntensity = pair.peaks[i, j].GetDenormalizedIntensity(pair.injectionTime);
                                            completeTotalIntensity[i, j] += denormalizedIntensity;
                                            completeTotalIntensity[i, numIsotopes] += denormalizedIntensity;
                                        }
                                        final[c, 1]++;
                                    }
                                }
                            } // End cluster loop
                        } // End isotope loop
                    } // End pair loop
                } // End pair list loop
            }
            if (allPairs != null && allPairs.Count > 0)
            {
                foreach (List<Pair> pairs in allPairs.Values)
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
                                            double denormalizedIntensity = pair.peaks[i, j].GetDenormalizedIntensity(pair.injectionTime);
                                            totalIntensity[i, j] += denormalizedIntensity;
                                            totalIntensity[i, numIsotopes] += denormalizedIntensity;

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

                                    if (noNullPeaks)
                                    {
                                        for (int i = channelIndex; i < channelIndex + numIsotopologues; i++)
                                        {
                                            if (pair.peaks[i, j] != null)
                                            {
                                                double denormalizedIntensity = pair.peaks[i, j].GetDenormalizedIntensity(pair.injectionTime);
                                                totalIntensity[i, j] += denormalizedIntensity;
                                                totalIntensity[i, numIsotopes] += denormalizedIntensity;
                                            }
                                        }
                                        final[c, 0]++;
                                    }
                                }
                            }
                        }
                    }
                }
            } 

            double lightInt;
            double heavyInt;
            bool quantified = false;
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
                            lightInt = completeTotalIntensity[index1, numIsotopes];
                            heavyInt = completeTotalIntensity[index2, numIsotopes];
                            finalQuantified[c] = final[c, 1];
                            if (lightInt > 0 && heavyInt > 0)
                            {
                                heavyToLightRatioSum[n - 1, 0] = heavyInt / lightInt;
                                quantifiedNoiseIncluded[c] = false;
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
                    int channelIndex = c * numIsotopologues;

                    if (final[c, 1] >= minimumPostQFPairs) //&& ((double)countNoNBCIsotopes[c] / (double)countAllIsotopes[c]) > 0.25)
                    {
                        finalQuantified[c] = final[c, 1];
                        // Use only complete pairs for quantification
                        for (int i = channelIndex ; i < channelIndex + numIsotopologues; i++)
                        {
                            int index1 = i;
                            int index2 = i + 1;
                            // Use only complete pairs for quantification
                            for (int n = index2; n < channelIndex + numIsotopologues; n++)
                            {
                                lightInt = completeTotalIntensity[index1, numIsotopes];
                                heavyInt = completeTotalIntensity[index2, numIsotopes];
                                finalQuantified[c] = final[c, 1];
                                if (lightInt > 0 && heavyInt > 0)
                                {
                                    heavyToLightRatioSum[n - 1, 0] = heavyInt / lightInt;
                                    quantifiedNoiseIncluded[c] = false;
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
                                quantifiedNoiseIncluded[c] = true;
                            }
                            //else if (Form1.NEUCODE_ARG_PROLINECONVERSION && lightInt > 0 && heavyInt == 0)
                            //{
                            //    heavyToLightRatioSum[i - 1, 0] = 0;
                            //    quantifiedNoiseIncluded[c] = true;
                            //}
                            //index2++;
                        }
                    }
                    else
                    {
                        // Not able to quantify
                        for (int i = channelIndex; i < channelIndex + numIsotopologues; i++)
                        {
                            //heavyToLightRatioSum[i - 1, 0] = double.NaN;
                            completeTotalIntensity[i, numIsotopes] = 0;
                            totalIntensity[i, numIsotopes] = 0;
                        }
                    }
                }
            }
        }
    }
}
