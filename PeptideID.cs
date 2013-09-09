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
        //public static MassTolerance SILAC = new MassTolerance(MassToleranceType.PPM, 20.0); // Individual SILAC peak tolerance (for lower resolution MS1 scans)
        //public static MassTolerance NEUCODE = new MassTolerance(MassToleranceType.PPM, 10.0); // Individual NeuCode peak tolerance (for higher resolution MS1 scans)
        public static double INTENSITYCUTOFF = 1.0 / (2.0 * Math.E); // Intensity threshold for quantitation filtering to eliminate low-level peaks
        public static double SIGNALTONOISE = Form1.MINIMUMSN;
        public static double QUANTRESOLUTION = Form1.QUANTRESOLUTION * 1000.0;
        public static double QUANTSEPARATION = Form1.THEORETICALSEPARATION;
        public static double SPACINGPERCENTAGEERROR = 0.50;
        public static MassTolerance TOLERANCE = new MassTolerance(MassToleranceType.PPM, (Form1.TOLERANCE * 2.0)); // User-specified ppm tolerance

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
        public int numIsotopologueLabels; // The number of quantitative labels carried by the peptide
        public int numClusterLabels;
        public int numLabels
        {
            get
            {
                if (numIsotopologues > 1) return numIsotopologueLabels;
                else return numClusterLabels;
            }
        }
        public string sequence; // A peptide's sequence
        public string sequenceNoMods; 
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
                int channels = numChannels;
                //if (Form1.NEUCODE_SIXPLEX_ARG)
                //{
                //    channels = numChannels + 2;
                //}
                //else
                //{
                //    channels = numChannels;
                //}

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
        public double[] averageTheoMZ
        {
            get
            {
                double[] averageTheoMZ = null;
                if (bestPSM != null)
                {
                    int charge = bestPSM.Charge;
                    if (numIsotopologues > 1)
                    {
                        averageTheoMZ = new double[numClusters];
                        for (int c = 0; c < numClusters; c++)
                        {
                            int channelIndex = c * numIsotopologues;

                            double avgTheoMass;
                            if (numIsotopologues % 2 == 0)
                            {
                                avgTheoMass = (theoMasses[channelIndex + ((numIsotopologues / 2) - 1), 0] + theoMasses[channelIndex + (numIsotopologues / 2), 0]) / 2.0;
                            }
                            else
                            {
                                avgTheoMass = theoMasses[channelIndex + (numIsotopologues / 2), 0];
                            }
                            averageTheoMZ[c] = Mass.MzFromMass(avgTheoMass, charge);
                        }
                        return averageTheoMZ;
                    }
                    else
                    {
                        averageTheoMZ = new double[numChannels];
                        for (int i = 0; i < numChannels; i++)
                        {
                            averageTheoMZ[i] = Mass.MzFromMass(theoMasses[i, 0], charge);
                        }
                        return averageTheoMZ;
                    }
                }
                return averageTheoMZ;
            }
        }
        public double[] averageAdjustedTheoMZ
        {
            get
            {
                double[] averageAdjustedTheoMZ = null;
                if (bestPSM != null)
                {
                    int charge = bestPSM.Charge;
                    if (numIsotopologues > 1)
                    {
                        averageAdjustedTheoMZ = new double[numClusters];
                        for (int c = 0; c < numClusters; c++)
                        {
                            int channelIndex = c * numIsotopologues;

                            double avgAdjTheoMass;
                            if (numIsotopologues % 2 == 0)
                            {
                                avgAdjTheoMass = (adjustedTheoMasses[channelIndex + ((numIsotopologues / 2) - 1), 0] + adjustedTheoMasses[channelIndex + (numIsotopologues / 2), 0]) / 2.0;
                            }
                            else
                            {
                                avgAdjTheoMass = adjustedTheoMasses[channelIndex + (numIsotopologues / 2), 0];
                            }
                            averageAdjustedTheoMZ[c] = Mass.MzFromMass(avgAdjTheoMass, charge);
                        }
                        return averageAdjustedTheoMZ;
                    }
                    else
                    {
                        averageAdjustedTheoMZ = new double[numChannels];
                        for (int i = 0; i < numChannels; i++)
                        {
                            averageAdjustedTheoMZ[i] = Mass.MzFromMass(adjustedTheoMasses[i, 0], charge);
                        }
                        return averageAdjustedTheoMZ;
                    }
                }
                return averageAdjustedTheoMZ;
            }
        }
        
        // Quantification information
        public NonQuantifiableType noQuantReason { get; set; }
        public int peaksNeeded
        {
            get
            {
                int numPeaks;
                if (Form1.NOISEBANDCAP)
                {
                    if (numIsotopologues < 2)
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
                    if (numIsotopologues < 2)
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
        //public MassTolerance tolerance
        //{
        //    get
        //    {
        //        MassTolerance tolerance;
        //        if (Form1.NEUCODE || Form1.SILAC_DUPLEX_LEUCN || Form1.SILAC_DUPLEX_LEUH)
        //        {
        //            tolerance = NEUCODE;
        //        }
        //        else
        //        {
        //            tolerance = SILAC;
        //        }
        //        return tolerance;
        //    }
        //}
        public double[,] maxIntensity // Considering all pairs, finds the intensity and retention time at which each peptide isotopologue reaches its elution maximum
        {
            get
            {
                double[,] max = null;
                List<Pair> pairsCombined;
                List<Pair> intensitySortedPairs;

                if (allPairs != null && allPairs.Count > 0)
                {
                    max = new double[numChannels, 2];
                    pairsCombined = new List<Pair>();
                    intensitySortedPairs = new List<Pair>();

                    foreach (List<Pair> pairList in allPairs.Values)
                    {
                        foreach (Pair pair in pairList)
                        {
                            pairsCombined.Add(pair);
                        }
                    }

                    for (int i = 0; i < numChannels; i++)
                    {
                        intensitySortedPairs = pairsCombined.OrderByDescending(channelPair => channelPair.GetMaxChannelIntensity(i)).Take(5).ToList();
                        switch (intensitySortedPairs.Count)
                        {
                            case 0:
                                max[i, 0] = 0;
                                max[i, 1] = 0;
                                break;
                            case 1:
                                max[i, 0] = intensitySortedPairs[0].GetMaxChannelIntensity(i);
                                max[i, 1] = intensitySortedPairs[0].retentionTime;
                                break;
                            case 2:
                                max[i, 0] = (intensitySortedPairs[0].GetMaxChannelIntensity(i) + intensitySortedPairs[1].GetMaxChannelIntensity(i)) / 2.0;
                                max[i, 1] = intensitySortedPairs[0].retentionTime;
                                break;
                            case 3:
                                max[i, 0] = intensitySortedPairs[1].GetMaxChannelIntensity(i);
                                max[i, 1] = intensitySortedPairs[0].retentionTime;
                                break;
                            case 4:
                                max[i, 0] = (intensitySortedPairs[1].GetMaxChannelIntensity(i) + intensitySortedPairs[2].GetMaxChannelIntensity(i)) / 2.0;
                                max[i, 1] = intensitySortedPairs[0].retentionTime;
                                break;
                            case 5:
                                max[i, 0] = intensitySortedPairs[2].GetMaxChannelIntensity(i);
                                max[i, 1] = intensitySortedPairs[0].retentionTime;
                                break;
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
                List<Pair> pairsCombined;
                List<Pair> intensitySortedPairs;

                if (completePairs != null && completePairs.Count > 0)
                {
                    max = new double[numChannels, 2];
                    pairsCombined = new List<Pair>();
                    intensitySortedPairs = new List<Pair>();

                    foreach (List<Pair> pairList in completePairs.Values)
                    {
                        foreach (Pair pair in pairList)
                        {
                            pairsCombined.Add(pair);
                        }
                    }

                    for (int i = 0; i < numChannels; i++)
                    {
                        intensitySortedPairs = pairsCombined.OrderByDescending(channelPair => channelPair.GetMaxChannelIntensity(i)).Take(5).ToList();
                        switch (intensitySortedPairs.Count)
                        {
                            case 0:
                                max[i, 0] = 0;
                                max[i, 1] = 0;
                                break;
                            case 1:
                                max[i, 0] = intensitySortedPairs[0].GetMaxChannelIntensity(i);
                                max[i, 1] = intensitySortedPairs[0].retentionTime;
                                break;
                            case 2:
                                max[i, 0] = (intensitySortedPairs[0].GetMaxChannelIntensity(i) + intensitySortedPairs[1].GetMaxChannelIntensity(i)) / 2.0;
                                max[i, 1] = intensitySortedPairs[0].retentionTime;
                                break;
                            case 3:
                                max[i, 0] = intensitySortedPairs[1].GetMaxChannelIntensity(i);
                                max[i, 1] = intensitySortedPairs[0].retentionTime;
                                break;
                            case 4:
                                max[i, 0] = (intensitySortedPairs[1].GetMaxChannelIntensity(i) + intensitySortedPairs[2].GetMaxChannelIntensity(i)) / 2.0;
                                max[i, 1] = intensitySortedPairs[0].retentionTime;
                                break;
                            case 5:
                                max[i, 0] = intensitySortedPairs[2].GetMaxChannelIntensity(i);
                                max[i, 1] = intensitySortedPairs[0].retentionTime;
                                break;
                        }
                    }
                }
                return max;
            }
        }
        public double[] log10MaxIntensity // Considering all isotopologues from all pairs, finds the peptide's maximum intensity
        {
            get
            {
                double[] max = null; ;
                if (absoluteMaxIntensity != null)
                {
                    if (numIsotopologues > 1)
                    {
                        max = new double[numClusters];
                        int channelIndex;
                        for (int c = 0; c < numClusters; c++)
                        {
                            channelIndex = c * numIsotopologues;
                            for (int i = channelIndex; i < channelIndex + numIsotopologues; i++)
                            {
                                if (absoluteMaxIntensity[c, 0] > max[c]) max[c] = absoluteMaxIntensity[c, 0];
                            }
                            max[c] = Math.Log10(max[c]);
                        }
                    }
                    else
                    {
                        max = new double[1];
                        for (int c = 0; c < numChannels; c++)
                        {
                            if (absoluteMaxIntensity[c, 0] > max[0]) max[0] = absoluteMaxIntensity[c, 0];
                        }
                        max[0] = Math.Log10(max[0]);
                    }
                }
                return max;
            }
        }
        public double[] log10MaxCompleteIntensity // Considering all isotopologues from complete pairs, finds the peptide's maximum intensity
        {
            get
            {
                double[] max = null; ;
                if (absoluteMaxCompleteIntensity != null)
                {
                    if (numIsotopologues > 1)
                    {
                        max = new double[numClusters];
                        int channelIndex;
                        for (int c = 0; c < numClusters; c++)
                        {
                            channelIndex = c * numIsotopologues;
                            for (int i = channelIndex; i < channelIndex + numIsotopologues; i++)
                            {
                                if (absoluteMaxCompleteIntensity[c, 0] > max[c]) max[c] = absoluteMaxCompleteIntensity[c, 0];
                            }
                            max[c] = Math.Log10(max[c]);
                        }
                    }
                    else
                    {
                        max = new double[1];
                        for (int c = 0; c < numChannels; c++)
                        {
                            if (absoluteMaxCompleteIntensity[c, 0] > max[0]) max[0] = absoluteMaxCompleteIntensity[c, 0];
                        }
                        max[0] = Math.Log10(max[0]);
                    }
                }
                return max;
            }
        }
        public double[] missingChannelFrequency
        {
            get
            {
                List<Pair> all = new List<Pair>();
                double[] frequency = new double[numClusters];
                int[] incompleteCount = new int[numClusters];
                int[] completeCount = new int[numClusters];
                int[] totalCount = new int[numClusters];
                for (int c = 0; c < numClusters; c++)
                {
                    frequency[c] = -1;
                    if (allPairs != null && allPairs.Count > 0) incompleteCount[c] = countAllIsotopes[c];
                    if (completePairs != null && completePairs.Count > 0) completeCount[c] = countCompleteIsotopes[c];
                    totalCount[c] = incompleteCount[c] + completeCount[c];

                    if (totalCount[c] > 0) frequency[c] = (double)incompleteCount[c] / (double)totalCount[c];
                }
                return frequency;
            }
        } // Calculate's a peptide's missing channel frequency to be used in coalescence detection
        public double[,] absoluteMaxIntensity
        {
            get
            {
                double[,] max = null;
                List<Pair> pairsCombined;
                List<Pair> intensitySortedPairs;

                if (allPairs != null && allPairs.Count > 0)
                {
                    max = new double[numChannels, 2];
                    pairsCombined = new List<Pair>();
                    intensitySortedPairs = new List<Pair>();

                    foreach (List<Pair> pairList in allPairs.Values)
                    {
                        foreach (Pair pair in pairList)
                        {
                            pairsCombined.Add(pair);
                        }
                    }

                    for (int i = 0; i < numChannels; i++)
                    {
                        intensitySortedPairs = pairsCombined.OrderByDescending(channelPair => channelPair.GetMaxChannelIntensity(i)).Take(1).ToList();
                        switch (intensitySortedPairs.Count)
                        {
                            case 0:
                                max[i, 0] = 0;
                                max[i, 1] = 0;
                                break;
                            case 1:
                                max[i, 0] = intensitySortedPairs[0].GetMaxChannelIntensity(i);
                                max[i, 1] = intensitySortedPairs[0].retentionTime;
                                break;
                        }
                    }
                }
                return max;
            }
        }
        public double[,] absoluteMaxCompleteIntensity
        {
            get
            {
                double[,] max = null;
                List<Pair> pairsCombined;
                List<Pair> intensitySortedPairs;

                if (completePairs != null && completePairs.Count > 0)
                {
                    max = new double[numChannels, 2];
                    pairsCombined = new List<Pair>();
                    intensitySortedPairs = new List<Pair>();

                    foreach (List<Pair> pairList in completePairs.Values)
                    {
                        foreach (Pair pair in pairList)
                        {
                            pairsCombined.Add(pair);
                        }
                    }

                    for (int i = 0; i < numChannels; i++)
                    {
                        intensitySortedPairs = pairsCombined.OrderByDescending(channelPair => channelPair.GetMaxChannelIntensity(i)).Take(1).ToList();
                        switch (intensitySortedPairs.Count)
                        {
                            case 0:
                                max[i, 0] = 0;
                                max[i, 1] = 0;
                                break;
                            case 1:
                                max[i, 0] = intensitySortedPairs[0].GetMaxChannelIntensity(i);
                                max[i, 1] = intensitySortedPairs[0].retentionTime;
                                break;
                        }
                    }
                }
                return max;
            }
        }
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
                    if (numIsotopologues > 1)
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
                                        if (channelCount >= peaksNeeded)
                                        {
                                            count[c]++;
                                        }
                                    }
                                }
                            }
                        }
                    }
                    else
                    {
                        count = new int[1];
                        foreach (MSDataFile rawFile in allPairs.Keys)
                        {
                            List<Pair> pairs = allPairs[rawFile];
                            foreach (Pair pair in pairs)
                            {
                                for (int j = 0; j < numIsotopes; j++)
                                {
                                    channelCount = 0;
                                    for (int i = 0; i < numChannels; i++)
                                    {
                                        if (pair.peaks[i, j] != null)
                                        {
                                            channelCount++;
                                        }
                                    }
                                    if (channelCount >= peaksNeeded)
                                    {
                                        count[0]++;
                                    }
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
        public int[,] missingChannelPattern; 
        public int[] finalQuantified;
        public MassRange[,] spacingMassRange
        {
            get
            {
                MassRange[,] spacing;
                double theoSpacing;
                double modificationSpacing;

                if (numIsotopologues > 1 && numIsotopologueLabels > 0)
                {
                    spacing = new MassRange[numIsotopologues, numIsotopologues];
                    theoSpacing = 0;

                    for (int r = 0; r < numIsotopologues; r++)
                    {
                        for (int c = 0; c < numIsotopologues; c++)
                        {
                            if (r == c) spacing[r, c] = null;
                            else
                            {
                                modificationSpacing = Form1.ISOTOPOLOGUELABELS[r].Mass.Monoisotopic - Form1.ISOTOPOLOGUELABELS[c].Mass.Monoisotopic;
                                theoSpacing = (double)numIsotopologueLabels * modificationSpacing;
                                if (theoSpacing < 0) theoSpacing = theoSpacing * -1.0;
                                spacing[r, c] = new MassRange(theoSpacing, new MassTolerance(MassToleranceType.DA, SPACINGPERCENTAGEERROR * theoSpacing));
                            }
                        }     
                    }
                    return spacing;
                }
                else if (numIsotopologues < 2 && numClusterLabels > 0)
                {
                    spacing = new MassRange[numChannels, numChannels];
                    theoSpacing = 0;

                    for (int r = 0; r < numChannels; r++)
                    {
                        for (int c = 0; c < numChannels; c++)
                        {
                            if (r == c) spacing[r, c] = null;
                            else if (r == 0)
                            {
                                theoSpacing = 0 - ((double)numClusterLabels * Form1.CLUSTERLABELS[c - 1].Mass.Monoisotopic);
                                if (theoSpacing < 0) theoSpacing = theoSpacing * -1.0;
                                spacing[r, c] = new MassRange(theoSpacing, new MassTolerance(MassToleranceType.DA, 0.020 * numClusterLabels));
                            }
                            else if (c == 0)
                            {
                                theoSpacing = ((double)numClusterLabels * Form1.CLUSTERLABELS[r - 1].Mass.Monoisotopic) - 0;
                                if (theoSpacing < 0) theoSpacing = theoSpacing * -1.0;
                                spacing[r, c] = new MassRange(theoSpacing, new MassTolerance(MassToleranceType.DA, 0.020 * numClusterLabels));
                            }
                            else
                            {
                                theoSpacing = (double)numClusterLabels * (Form1.CLUSTERLABELS[r - 1].Mass.Monoisotopic - Form1.CLUSTERLABELS[c - 1].Mass.Monoisotopic);
                                if (theoSpacing < 0) theoSpacing = theoSpacing * -1.0;
                                spacing[r, c] = new MassRange(theoSpacing, new MassTolerance(MassToleranceType.DA, 0.020 * numClusterLabels));
                            }
                        }
                    }
                    return spacing;
                }
                else
                {
                    return null;
                }


                //if (numLabels > 0)
                //{
                //    if (numIsotopologues > 1)
                //    {
                //        spacing = new MassRange[numIsotopologues - 1];
                //        theoSpacing = new double[numIsotopologues - 1];

                //        for (int i = 1; i <= numIsotopologues; i++)
                //        {
                //            modificationSpacing = Form1.ISOTOPOLOGUELABELS[i].Mass.Monoisotopic - Form1.ISOTOPOLOGUELABELS[i-1].Mass.Monoisotopic;
                //            theoSpacing[i - 1] = (double)numLabels * modificationSpacing;
                //            spacing[i - 1] = new MassRange(theoSpacing[i - 1], new MassTolerance(MassToleranceType.DA, SPACINGPERCENTAGEERROR * theoSpacing[i - 1]));
                //        }                        
                        
                        //if (numIsotopologues == 6)
                        //{
                        //    theoSpacing1 = (double)numLabels * (2 * (Constants.DEUTERIUM - Constants.HYDROGEN) - 2 * (Constants.CARBON13 - Constants.CARBON));
                        //    theoSpacing2 = (double)numLabels * ((Constants.CARBON13 - Constants.CARBON) - (Constants.NITROGEN15 - Constants.NITROGEN));
                        //    theoSpacing3 = (double)numLabels * (4 * (Constants.DEUTERIUM - Constants.HYDROGEN) + (Constants.NITROGEN15 - Constants.NITROGEN) - 5 * (Constants.CARBON13 - Constants.CARBON));
                        //    theoSpacing4 = (double)numLabels * (4 * (Constants.CARBON13 - Constants.CARBON) - 2 * (Constants.DEUTERIUM - Constants.HYDROGEN) - 2 * (Constants.NITROGEN15 - Constants.NITROGEN));
                        //    theoSpacing5 = (double)numLabels * (4 * (Constants.DEUTERIUM - Constants.HYDROGEN) - 4 * (Constants.CARBON13 - Constants.CARBON));
                        //    spacing[0] = new MassRange(theoSpacing1, new MassTolerance(MassToleranceType.DA, SPACINGPERCENTAGEERROR * theoSpacing1));
                        //    spacing[1] = new MassRange(theoSpacing2, new MassTolerance(MassToleranceType.DA, SPACINGPERCENTAGEERROR * theoSpacing2));
                        //    spacing[2] = new MassRange(theoSpacing3, new MassTolerance(MassToleranceType.DA, SPACINGPERCENTAGEERROR * theoSpacing3));
                        //    spacing[3] = new MassRange(theoSpacing4, new MassTolerance(MassToleranceType.DA, SPACINGPERCENTAGEERROR * theoSpacing4));
                        //    spacing[4] = new MassRange(theoSpacing5, new MassTolerance(MassToleranceType.DA, SPACINGPERCENTAGEERROR * theoSpacing5));
                        //}
                        //else if (numIsotopologues == 4)
                        //{
                        //    if (Form1.NEUCODE_FOURPLEX_LYS8_12MDA)
                        //    {
                        //        theoSpacing1 = (double)numLabels * (2 * (Constants.DEUTERIUM - Constants.HYDROGEN) - (Constants.NITROGEN15 - Constants.NITROGEN) - (Constants.CARBON13 - Constants.CARBON));
                        //        theoSpacing2 = theoSpacing1;
                        //        theoSpacing3 = (double)numLabels * (4 * (Constants.DEUTERIUM - Constants.HYDROGEN) - 4 * (Constants.CARBON13 - Constants.CARBON));
                        //        spacing[0] = new MassRange(theoSpacing1, new MassTolerance(MassToleranceType.DA, SPACINGPERCENTAGEERROR * theoSpacing1));
                        //        spacing[1] = spacing[0];
                        //        spacing[2] = new MassRange(theoSpacing3, new MassTolerance(MassToleranceType.DA, SPACINGPERCENTAGEERROR * theoSpacing3));
                        //    }
                        //    else
                        //    {
                        //        theoSpacing1 = (double)numLabels * (2 * (Constants.CARBON13 - Constants.CARBON) - 2 * (Constants.NITROGEN15 - Constants.NITROGEN));
                        //        theoSpacing2 = theoSpacing1;
                        //        theoSpacing3 = theoSpacing1;
                        //        spacing[0] = new MassRange(theoSpacing1, new MassTolerance(MassToleranceType.DA, SPACINGPERCENTAGEERROR * theoSpacing1));
                        //        spacing[1] = spacing[0];
                        //        spacing[2] = spacing[1];
                        //    }
                        //}
                        //else if (numIsotopologues == 3)
                        //{
                        //    theoSpacing1 = (double)numLabels * (6 * (Constants.DEUTERIUM - Constants.HYDROGEN) - 6 * (Constants.CARBON13 - Constants.CARBON));
                        //    theoSpacing2 = (double)numLabels * (2 * (Constants.DEUTERIUM - Constants.HYDROGEN) - 2 * (Constants.NITROGEN15 - Constants.NITROGEN));
                        //    spacing[0] = new MassRange(theoSpacing1, new MassTolerance(MassToleranceType.DA, SPACINGPERCENTAGEERROR * theoSpacing1));
                        //    spacing[1] = new MassRange(theoSpacing2, new MassTolerance(MassToleranceType.DA, SPACINGPERCENTAGEERROR * theoSpacing2));
                        //}
                        //else if (Form1.NEUCODE_DUPLEX_LYS1_6MDA || Form1.NEUCODE_DUPLEX_CARBAMYL)
                        //{
                        //    theoSpacing1 = (double)numLabels * ((Constants.CARBON13 - Constants.CARBON) - (Constants.NITROGEN15 - Constants.NITROGEN));
                        //    spacing[0] = new MassRange(theoSpacing1, new MassTolerance(MassToleranceType.DA, SPACINGPERCENTAGEERROR * theoSpacing1));
                        //}
                        //else if (Form1.NEUCODE_DUPLEX_LEU7_18MDA)
                        //{
                        //    if (conversionFactor > 0)
                        //    {
                        //        theoSpacing1 = ((double)numLabels * (7 * (Constants.DEUTERIUM - Constants.HYDROGEN) - (6 * (Constants.CARBON13 - Constants.CARBON) + (1 * (Constants.NITROGEN15 - Constants.NITROGEN))))) - conversionFactor * (Constants.NITROGEN15 - Constants.NITROGEN);
                        //    }
                        //    else
                        //    {
                        //        theoSpacing1 = (double)numLabels * (7 * (Constants.DEUTERIUM - Constants.HYDROGEN) - (6 * (Constants.CARBON13 - Constants.CARBON) + (1 * (Constants.NITROGEN15 - Constants.NITROGEN))));
                        //    }
                        //    spacing[0] = new MassRange(theoSpacing1, new MassTolerance(MassToleranceType.DA, SPACINGPERCENTAGEERROR * theoSpacing1));
                        //}
                        //else
                        //{
                        //    theoSpacing1 = (double)numLabels * (8 * (Constants.DEUTERIUM - Constants.HYDROGEN) - (6 * (Constants.CARBON13 - Constants.CARBON) + (2 * (Constants.NITROGEN15 - Constants.NITROGEN))));
                        //    spacing[0] = new MassRange(theoSpacing1, new MassTolerance(MassToleranceType.DA, SPACINGPERCENTAGEERROR * theoSpacing1));
                        //}

                    //}
                    //else
                    //{
                        //spacing = new MassRange[numClusters - 1];
                        //theoSpacing = new double[numClusters - 1];

                        //theoSpacing[0] = (double)numLabels * Form1.CLUSTERLABELS[0].Mass.Monoisotopic;
                        //spacing[0] = new MassRange(theoSpacing[0], new MassTolerance(MassToleranceType.DA, 0.010 * numLabels));
                        
                        //if (numClusters == 3)
                        //{
                        //    theoSpacing[1] = (double)numLabels * (Form1.CLUSTERLABELS[1].Mass.Monoisotopic - Form1.CLUSTERLABELS[0].Mass.Monoisotopic);
                        //    spacing[1] = new MassRange(theoSpacing[1], new MassTolerance(MassToleranceType.DA, 0.01 * numLabels));
                        //}

                        //if (Form1.SILAC_DUPLEX_LYSC)
                        //{
                        //    theoSpacing1 = numLabels * (6 * (Constants.CARBON13 - Constants.CARBON));
                        //}                        
                        //else if (Form1.SILAC_DUPLEX_LYSH)
                        //{
                        //    theoSpacing1 = numLabels * (8 * (Constants.DEUTERIUM - Constants.HYDROGEN));
                        //}
                        //else if (Form1.SILAC_DUPLEX_LEUCN)
                        //{
                        //    if (conversionFactor > 0)
                        //    {
                        //        theoSpacing1 = numLabels * ((6 * (Constants.CARBON13 - Constants.CARBON)) + (Constants.NITROGEN15 - Constants.NITROGEN)) - conversionFactor * (Constants.NITROGEN15 - Constants.NITROGEN);
                        //    }
                        //    else
                        //    {
                        //        theoSpacing1 = numLabels * ((6 * (Constants.CARBON13 - Constants.CARBON)) + (Constants.NITROGEN15 - Constants.NITROGEN));
                        //    }
                        //}
                        //else if (Form1.SILAC_DUPLEX_LEUH)
                        //{
                        //    theoSpacing1 = numLabels * (7 * (Constants.DEUTERIUM - Constants.HYDROGEN));
                        //}
                        //else
                        //{
                        //    theoSpacing1 = numLabels * ((6 * (Constants.CARBON13 - Constants.CARBON)) + (2 * (Constants.NITROGEN15 - Constants.NITROGEN)));
                        //}
                        //spacing[0] = new MassRange(theoSpacing1, new MassTolerance(MassToleranceType.DA, 0.010 * numLabels));
                    //}
                //    return spacing;
                //}
                //else
                //{
                //    return null;
                //}
            }
        }
        public int conversionFactor { get; set; }
        public bool theoreticallyResolvable
        {
            get
            {
                if (numIsotopologues < 2)
                {
                    if (numLabels > 0) return true;
                    else return false;
                }
                else
                {
                    for (int i = 0; i < numClusters; i++)
                    {
                        if (GetTheoreticalResolvability(i))
                        {
                            return true;
                        }
                    }
                    return false;
                }
            }
        }

        // Coalescence information
        public bool coalescenceDetected;
        public List<double> coalescedPeakIntensities;
        public Dictionary<int, List<double>> missingChannelsSN;

        // Creates a PeptideID based on a database search peptide identification
        public PeptideID(int scanNumber, int charge, double eValue, string sequenceOriginal, MSDataFile rawFile, string mods)
        {           
            // Local variables
            List<PeptideSpectralMatch> psms;
            Peptide peptide;
            Peptide[] peptideVersions = null;

            // Get rid of variable label incorporations
            string sequenceFixed = "";
            for (int i = 0; i < sequenceOriginal.Length; i++)
            {
                if (sequenceOriginal[i].Equals('l')) sequenceFixed += 'L';
                else if (sequenceOriginal[i].Equals('k') && !mods.Contains("q")) sequenceFixed += 'K';
                else sequenceFixed += sequenceOriginal[i];
            }
            
            // Initialize all peptide properties
            sequence = sequenceFixed;
            numChannels = Form1.NUMCHANNELS;
            numIsotopes = Form1.NUMISOTOPES;
            numIsotopologues = Form1.NUMISOTOPOLOGUES;
            numClusters = Form1.NUMCLUSTERS;

            // Deal with variable mods that do not affect quantification (e.g., oxidation, phosphorylation)
            sequenceNoMods = "";
            List<int> oxidationPositions = new List<int>();
            List<int> phosphorylationPositions = new List<int>();
            List<int> tyrosineNHSPositions = new List<int>();
            List<int> kGGPositions = new List<int>();
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
                    if (Form1.NHSCLUSTER)
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
                    kGGPositions.Add(i);
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
                    peptide.SetModification(NamedChemicalFormula.Phosphorylation, position + 1);
                }
            }

            //missingChannelsSN = new Dictionary<int, List<double>>();
            //missingChannelsSN.Add(1, new List<double>());
            //missingChannelsSN.Add(2, new List<double>());
            //missingChannelsSN.Add(3, new List<double>());
            //missingChannelsSN.Add(4, new List<double>());

            // Set number of labels and search range for PPM correction
            if (numIsotopologues > 1)
            {
                if (Form1.NHSISOTOPOLOGUE)
                {
                    numIsotopologueLabels = countResidues('K', peptide.Sequence) + 1 + tyrosineNHSPositions.Count;
                    firstSearchMassRange = new MassTolerance(MassToleranceType.PPM, 50.0);
                }
                else if (Form1.LEUISOTOPOLOGUE)
                {
                    numIsotopologueLabels = countResidues('L', peptide.Sequence);
                    firstSearchMassRange = new MassTolerance(MassToleranceType.PPM, 50.0);
                }
                else
                {
                    numIsotopologueLabels = countResidues('K', peptide.Sequence);
                    firstSearchMassRange = new MassTolerance(MassToleranceType.PPM, 50.0);

                    if (numClusters > 1)
                    {
                        if (Form1.NHSCLUSTER)
                        {
                            numClusterLabels = countResidues('K', peptide.Sequence) + 1 + tyrosineNHSPositions.Count;
                        }
                    }
                }
            }
            else
            {
                if (Form1.LEUCLUSTER)
                {
                    numClusterLabels = countResidues('L', peptide.Sequence);
                    firstSearchMassRange = new MassTolerance(MassToleranceType.PPM, 50.0);
                }
                else if (Form1.CYSCLUSTER)
                {
                    numClusterLabels = countResidues('C', peptide.Sequence);
                    firstSearchMassRange = new MassTolerance(MassToleranceType.PPM, 50.0);
                }
                else if (Form1.ARGCLUSTER)
                {
                    numClusterLabels = countResidues('R', peptide.Sequence);
                    firstSearchMassRange = new MassTolerance(MassToleranceType.PPM, 50.0);
                }
                else if (Form1.NHSCLUSTER)
                {
                    numClusterLabels = countResidues('K', peptide.Sequence) + 1 + tyrosineNHSPositions.Count;
                    firstSearchMassRange = new MassTolerance(MassToleranceType.PPM, 50.0);
                }
                else
                {
                    numClusterLabels = countResidues('K', peptide.Sequence);
                    firstSearchMassRange = new MassTolerance(MassToleranceType.PPM, 50.0);
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
            if (numIsotopologues > 1)
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
            missingChannelPattern = new int[numClusters, 2];
            int numPeptideVersions;
            
            // Set theoretical masses of each channel
            if (numIsotopologues > 1 && numIsotopologueLabels > 0)
            {
                if (numClusters > 1)
                {
                    if (numClusterLabels > 0) numPeptideVersions = numIsotopologues * numClusters;
                    else numPeptideVersions = numIsotopologues;
                }
                else
                {
                    numPeptideVersions = numIsotopologues;
                }
                peptideVersions = new Peptide[numPeptideVersions];
                
                for (int i = 0; i < numPeptideVersions; i++)
                {
                    peptideVersions[i] = new Peptide(peptide);
                }

                if ((numClusters == 1 && Form1.CLUSTERLABELS.Count == 0) || numClusterLabels == 0)
                {
                    for (int i = 0; i < Form1.ISOTOPOLOGUELABELS.Count; i++)
                    {
                        if (Form1.LYSISOTOPOLOGUE)
                        {
                            peptideVersions[i].SetModification(Form1.ISOTOPOLOGUELABELS[i], ModificationSites.K);

                            if (kGGPositions.Count > 0)
                            {
                                foreach (int position in kGGPositions)
                                {
                                    peptideVersions[i].SetModification(new Mass(Form1.ISOTOPOLOGUELABELS[i].Mass.Monoisotopic + 114.0429), position + 1);
                                }
                            }
                        }
                        else if (Form1.LEUISOTOPOLOGUE)
                        {
                            peptideVersions[i].SetModification(Form1.ISOTOPOLOGUELABELS[i], ModificationSites.L);

                            if (kGGPositions.Count > 0)
                            {
                                foreach (int position in kGGPositions)
                                {
                                    peptideVersions[i].SetModification(new Mass(114.0429), position + 1);
                                }
                            }
                        }
                        else if (Form1.NHSISOTOPOLOGUE)
                        {
                            peptideVersions[i].SetModification(Form1.ISOTOPOLOGUELABELS[i], ModificationSites.NPep | ModificationSites.K);

                            if (tyrosineNHSPositions.Count > 0)
                            {
                                foreach (int position in tyrosineNHSPositions)
                                {
                                    peptideVersions[i].SetModification(Form1.ISOTOPOLOGUELABELS[i], position);
                                }
                            }
                        }
                    }
                }

                else if (numClusters > 2)
                {
                    int channelIndex;

                    for (int c = 0; c < numClusters; c++)
                    {
                        channelIndex = c * numIsotopologues;
                        for (int i = channelIndex; i < channelIndex + numIsotopologues; i++)
                        {
                            if (Form1.ARGCLUSTER)
                            {
                                if (c > 0)
                                {
                                    peptideVersions[i].SetModification(Form1.CLUSTERLABELS[c - 1], ModificationSites.R);
                                }
                                    
                                if (Form1.LYSISOTOPOLOGUE)
                                {
                                    peptideVersions[i].SetModification(Form1.ISOTOPOLOGUELABELS[i], ModificationSites.K);

                                    if (kGGPositions.Count > 0)
                                    {
                                        foreach (int position in kGGPositions)
                                        {
                                            peptideVersions[i].SetModification(new Mass(Form1.ISOTOPOLOGUELABELS[i].Mass.Monoisotopic + 114.0429), position + 1);
                                        }
                                    }
                                }
                                else if (Form1.LEUISOTOPOLOGUE)
                                {
                                    peptideVersions[i].SetModification(Form1.ISOTOPOLOGUELABELS[i], ModificationSites.L);

                                    if (kGGPositions.Count > 0)
                                    {
                                        foreach (int position in kGGPositions)
                                        {
                                            peptideVersions[i].SetModification(new Mass(114.0429), position + 1);
                                        }
                                    }
                                }
                            }
                            else if (Form1.LEUCLUSTER)
                            {
                                if (c > 0)
                                {
                                    peptideVersions[i].SetModification(Form1.CLUSTERLABELS[c - 1], ModificationSites.L);
                                }
                                    
                                if (Form1.LYSISOTOPOLOGUE)
                                {
                                    peptideVersions[i].SetModification(Form1.ISOTOPOLOGUELABELS[i], ModificationSites.K);

                                    if (kGGPositions.Count > 0)
                                    {
                                        foreach (int position in kGGPositions)
                                        {
                                            peptideVersions[i].SetModification(new Mass(Form1.ISOTOPOLOGUELABELS[i].Mass.Monoisotopic + 114.0429), position + 1);
                                        }
                                    }
                                }
                            }
                            else if (Form1.NHSCLUSTER)
                            {
                                peptideVersions[i].SetModification(Form1.CLUSTERLABELS[c], ModificationSites.NPep);
                                peptideVersions[i].SetModification(Form1.ISOTOPOLOGUELABELS[i], ModificationSites.K);

                                if (tyrosineNHSPositions.Count > 0)
                                {
                                    foreach (int position in tyrosineNHSPositions)
                                    {
                                        peptideVersions[i].SetModification(Form1.CLUSTERLABELS[c], position);
                                    }
                                }
                            }
                        }
                    }
                }
            }
            else if (numIsotopologues < 2 && numClusterLabels > 0)
            {
                numPeptideVersions = numClusters;
                peptideVersions = new Peptide[numPeptideVersions];

                for (int i = 0; i < numPeptideVersions; i++)
                {
                    peptideVersions[i] = new Peptide(peptide);
                }

                for (int i = 0; i < Form1.CLUSTERLABELS.Count; i++)
                {
                    if (Form1.LYSCLUSTER)
                    {
                        peptideVersions[i + 1].SetModification(Form1.CLUSTERLABELS[i], ModificationSites.K);
                    }
                    else if (Form1.ARGCLUSTER)
                    {
                        peptideVersions[i + 1].SetModification(Form1.CLUSTERLABELS[i], ModificationSites.R);
                    }
                    else if (Form1.LEUCLUSTER)
                    {
                        peptideVersions[i + 1].SetModification(Form1.CLUSTERLABELS[i], ModificationSites.L);
                    }
                    else if (Form1.CYSCLUSTER)
                    {
                        peptideVersions[i + 1].SetModification(Form1.CLUSTERLABELS[i], ModificationSites.L);
                    }
                    else if (Form1.NHSCLUSTER)
                    {
                        peptideVersions[i].SetModification(Form1.CLUSTERLABELS[i], ModificationSites.NPep | ModificationSites.K);
                        
                        if (tyrosineNHSPositions.Count > 0)
                        {
                            foreach (int position in tyrosineNHSPositions)
                            {
                                peptideVersions[i].SetModification(Form1.CLUSTERLABELS[i], position);
                            }
                        }
                    
                    }
                }
            }

            if (peptideVersions != null)
            {
                for (int i = 0; i < peptideVersions.Length; i++)
                {
                    theoMasses[i, 0] = peptideVersions[i].Mass.Monoisotopic;
                }
            }
                
                // Duplex NeuCode (metabolic)
                //if (Form1.NEUCODE_DUPLEX_LYS8_36MDA || Form1.NEUCODE_DUPLEX_LYS1_6MDA || Form1.NEUCODE_DUPLEX_LEU7_18MDA)
                //{
                //    check1 = new Peptide(peptide);
                //    check2 = new Peptide(peptide);
                //    if (Form1.NEUCODE_DUPLEX_LYS8_36MDA)
                //    {
                //        check1.SetModification(NamedChemicalFormula.GetModification("Lys +8 13C6 15N2"), ModificationSites.K);
                //        check2.SetModification(NamedChemicalFormula.GetModification("Lys +8 2H8"), ModificationSites.K);
                //    }
                //    else if (Form1.NEUCODE_DUPLEX_LYS1_6MDA)
                //    {
                //        check1.SetModification(NamedChemicalFormula.GetModification("Lys +1 15N"), ModificationSites.K);
                //        check2.SetModification(NamedChemicalFormula.GetModification("Lys +1 13C"), ModificationSites.K);
                //    }
                //    theoMasses[0, 0] = check1.Mass.Monoisotopic;
                //    theoMasses[1, 0] = check2.Mass.Monoisotopic;
                //}

                // 3plex NeuCode
                //else if (Form1.NEUCODE_TRIPLEX_LYS8_18MDA)
                //{
                //    check1 = new Peptide(peptide);
                //    check2 = new Peptide(peptide);
                //    check3 = new Peptide(peptide);
                //    check1.SetModification(NamedChemicalFormula.GetModification("Lys +8 13C6 15N2"), ModificationSites.K);
                //    check2.SetModification(NamedChemicalFormula.GetModification("Lys +8 2H6 15N2"), ModificationSites.K);
                //    check3.SetModification(NamedChemicalFormula.GetModification("Lys +8 2H8"), ModificationSites.K);
                //    theoMasses[0, 0] = check1.Mass.Monoisotopic;
                //    theoMasses[1, 0] = check2.Mass.Monoisotopic;
                //    theoMasses[2, 0] = check3.Mass.Monoisotopic;
                //}

                // 4plex NeuCode (metabolic)
                //else if (Form1.NEUCODE_FOURPLEX_LYS8_12MDA)
                //{                    
                //    check1 = new Peptide(peptide);
                //    check2 = new Peptide(peptide);
                //    check3 = new Peptide(peptide);
                //    check4 = new Peptide(peptide);
                //    check1.SetModification(NamedChemicalFormula.GetModification("Lys +8 13C6 15N2"), ModificationSites.K);
                //    check2.SetModification(NamedChemicalFormula.GetModification("Lys +8 13C5 2H2 15N1"), ModificationSites.K);
                //    check3.SetModification(NamedChemicalFormula.GetModification("Lys +8 13C4 2H4"), ModificationSites.K);
                //    check4.SetModification(NamedChemicalFormula.GetModification("Lys +8 2H8"), ModificationSites.K);
                //    if (kGGPositions.Count > 0)
                //    {
                //        foreach (int position in kGGPositions)
                //        {
                //            check1.SetModification(new Mass(NamedChemicalFormula.GetModification("Lys +8 13C6 15N2").Mass.Monoisotopic + 114.0429), position + 1);
                //            check2.SetModification(new Mass(NamedChemicalFormula.GetModification("Lys +8 13C5 2H2 15N1").Mass.Monoisotopic + 114.0429), position + 1);
                //            check3.SetModification(new Mass(NamedChemicalFormula.GetModification("Lys +8 13C4 2H4").Mass.Monoisotopic + 114.0429), position + 1);
                //            check4.SetModification(new Mass(NamedChemicalFormula.GetModification("Lys +8 2H8").Mass.Monoisotopic + 114.0429), position + 1);
                //        }
                //    }

                //    theoMasses[0, 0] = check1.Mass.Monoisotopic;
                //    theoMasses[1, 0] = check2.Mass.Monoisotopic;
                //    theoMasses[2, 0] = check3.Mass.Monoisotopic;
                //    theoMasses[3, 0] = check4.Mass.Monoisotopic;
                //}

                // 6plex NeuCode (metabolic)
                //else if (Form1.NEUCODE_SIXPLEX_LYS8_6MDA)
                //{
                //    check1 = new Peptide(peptide);
                //    check2 = new Peptide(peptide);
                //    check3 = new Peptide(peptide);
                //    check4 = new Peptide(peptide);
                //    check5 = new Peptide(peptide);
                //    check6 = new Peptide(peptide);
                //    check1.SetModification(NamedChemicalFormula.GetModification("Lys +8 13C6 15N2"), ModificationSites.K);
                //    check2.SetModification(NamedChemicalFormula.GetModification("Lys +8 13C4 2H2 15N2"), ModificationSites.K);
                //    check3.SetModification(NamedChemicalFormula.GetModification("Lys +8 13C5 2H2 15N1"), ModificationSites.K);
                //    check4.SetModification(NamedChemicalFormula.GetModification("Lys +8 2H6 15N2"), ModificationSites.K);
                //    check5.SetModification(NamedChemicalFormula.GetModification("Lys +8 13C4 2H4"), ModificationSites.K);
                //    check6.SetModification(NamedChemicalFormula.GetModification("Lys +8 2H8"), ModificationSites.K);
                //    theoMasses[0, 0] = check1.Mass.Monoisotopic;
                //    theoMasses[1, 0] = check2.Mass.Monoisotopic;
                //    theoMasses[2, 0] = check3.Mass.Monoisotopic;
                //    theoMasses[3, 0] = check4.Mass.Monoisotopic;
                //    theoMasses[4, 0] = check5.Mass.Monoisotopic;
                //    theoMasses[5, 0] = check6.Mass.Monoisotopic;
                //}

                // Duplex NeuCode (chemical)
                //else if (Form1.NEUCODE_DUPLEX_CARBAMYL)
                //{
                //    check1 = new Peptide(peptide);
                //    check2 = new Peptide(peptide);

                //    check1.SetModification(NamedChemicalFormula.GetModification("Carbamyl L"), ModificationSites.K | ModificationSites.NPep);
                //    check2.SetModification(NamedChemicalFormula.GetModification("Carbamyl H"), ModificationSites.K | ModificationSites.NPep);
                //    if (tyrosineNHSPositions.Count > 0)
                //    {
                //        foreach (int position in tyrosineNHSPositions)
                //        {
                //            check1.SetModification(NamedChemicalFormula.GetModification("Carbamyl L"), position);
                //            check2.SetModification(NamedChemicalFormula.GetModification("Carbamyl H"), position);
                //        }
                //    }
                //    theoMasses[0, 0] = check1.Mass.Monoisotopic;
                //    theoMasses[1, 0] = check2.Mass.Monoisotopic;
                //}

                // 4plex NeuCode (chemical)
                //else if (Form1.NEUCODE_4PLEX_LIGHT || Form1.NEUCODE_4PLEX_MEDIUM || Form1.NEUCODE_4PLEX_HEAVY)
                //{
                //    check1 = new Peptide(peptide);
                //    check2 = new Peptide(peptide);
                //    check3 = new Peptide(peptide);
                //    check4 = new Peptide(peptide);
                //    string mod1;
                //    string mod2;
                //    string mod3;
                //    string mod4;
                //    if (Form1.NEUCODE_4PLEX_LIGHT)
                //    {
                //        mod1 = "4plex L1";
                //        mod2 = "4plex L2";
                //        mod3 = "4plex L3";
                //        mod4 = "4plex L4";
                //    }
                //    else if (Form1.NEUCODE_4PLEX_MEDIUM)
                //    {
                //        mod1 = "4plex M1";
                //        mod2 = "4plex M2";
                //        mod3 = "4plex M3";
                //        mod4 = "4plex M4";
                //    }
                //    else
                //    {
                //        mod1 = "4plex H1";
                //        mod2 = "4plex H2";
                //        mod3 = "4plex H3";
                //        mod4 = "4plex H4";
                //    }
                //    check1.SetModification(NamedChemicalFormula.GetModification(mod1), ModificationSites.K | ModificationSites.NPep);
                //    check2.SetModification(NamedChemicalFormula.GetModification(mod2), ModificationSites.K | ModificationSites.NPep);
                //    check3.SetModification(NamedChemicalFormula.GetModification(mod3), ModificationSites.K | ModificationSites.NPep);
                //    check4.SetModification(NamedChemicalFormula.GetModification(mod4), ModificationSites.K | ModificationSites.NPep);
                //    if (tyrosineNHSPositions.Count > 0)
                //    {
                //        foreach (int position in tyrosineNHSPositions)
                //        {
                //            check1.SetModification(NamedChemicalFormula.GetModification(mod1), position);
                //            check2.SetModification(NamedChemicalFormula.GetModification(mod2), position);
                //            check3.SetModification(NamedChemicalFormula.GetModification(mod3), position);
                //            check4.SetModification(NamedChemicalFormula.GetModification(mod4), position);
                //        }
                //    }
                //    theoMasses[0, 0] = check1.Mass.Monoisotopic;
                //    theoMasses[1, 0] = check2.Mass.Monoisotopic;
                //    theoMasses[2, 0] = check3.Mass.Monoisotopic;
                //    theoMasses[3, 0] = check4.Mass.Monoisotopic;
                //}

                // 12plex NeuCode (chemical)
                //else if (Form1.NEUCODE_12PLEX)
                //{
                //    check1 = new Peptide(peptide);
                //    check2 = new Peptide(peptide);
                //    check3 = new Peptide(peptide);
                //    check4 = new Peptide(peptide);
                //    check5 = new Peptide(peptide);
                //    check6 = new Peptide(peptide);
                //    check7 = new Peptide(peptide);
                //    check8 = new Peptide(peptide);
                //    check9 = new Peptide(peptide);
                //    check10 = new Peptide(peptide);
                //    check11 = new Peptide(peptide);
                //    check12 = new Peptide(peptide);

                //    check1.SetModification(NamedChemicalFormula.GetModification("4plex L1"), ModificationSites.K | ModificationSites.NPep);
                //    check2.SetModification(NamedChemicalFormula.GetModification("4plex L2"), ModificationSites.K | ModificationSites.NPep);
                //    check3.SetModification(NamedChemicalFormula.GetModification("4plex L3"), ModificationSites.K | ModificationSites.NPep);
                //    check4.SetModification(NamedChemicalFormula.GetModification("4plex L4"), ModificationSites.K | ModificationSites.NPep);
                //    check5.SetModification(NamedChemicalFormula.GetModification("4plex M1"), ModificationSites.K | ModificationSites.NPep);
                //    check6.SetModification(NamedChemicalFormula.GetModification("4plex M2"), ModificationSites.K | ModificationSites.NPep);
                //    check7.SetModification(NamedChemicalFormula.GetModification("4plex M3"), ModificationSites.K | ModificationSites.NPep);
                //    check8.SetModification(NamedChemicalFormula.GetModification("4plex M4"), ModificationSites.K | ModificationSites.NPep);
                //    check9.SetModification(NamedChemicalFormula.GetModification("4plex H1"), ModificationSites.K | ModificationSites.NPep);
                //    check10.SetModification(NamedChemicalFormula.GetModification("4plex H2"), ModificationSites.K | ModificationSites.NPep);
                //    check11.SetModification(NamedChemicalFormula.GetModification("4plex H3"), ModificationSites.K | ModificationSites.NPep);
                //    check12.SetModification(NamedChemicalFormula.GetModification("4plex H4"), ModificationSites.K | ModificationSites.NPep);
                //    if (tyrosineNHSPositions.Count > 0)
                //    {
                //        foreach (int position in tyrosineNHSPositions)
                //        {
                //            check1.SetModification(NamedChemicalFormula.GetModification("4plex L1"), position);
                //            check2.SetModification(NamedChemicalFormula.GetModification("4plex L2"), position);
                //            check3.SetModification(NamedChemicalFormula.GetModification("4plex L3"), position);
                //            check4.SetModification(NamedChemicalFormula.GetModification("4plex L4"), position);
                //            check5.SetModification(NamedChemicalFormula.GetModification("4plex M1"), position);
                //            check6.SetModification(NamedChemicalFormula.GetModification("4plex M2"), position);
                //            check7.SetModification(NamedChemicalFormula.GetModification("4plex M3"), position);
                //            check8.SetModification(NamedChemicalFormula.GetModification("4plex M4"), position);
                //            check9.SetModification(NamedChemicalFormula.GetModification("4plex H1"), position);
                //            check10.SetModification(NamedChemicalFormula.GetModification("4plex H2"), position);
                //            check11.SetModification(NamedChemicalFormula.GetModification("4plex H3"), position);
                //            check12.SetModification(NamedChemicalFormula.GetModification("4plex H4"), position);
                //        }
                //    }

                //    theoMasses[0, 0] = check1.Mass.Monoisotopic;
                //    theoMasses[1, 0] = check2.Mass.Monoisotopic;
                //    theoMasses[2, 0] = check3.Mass.Monoisotopic;
                //    theoMasses[3, 0] = check4.Mass.Monoisotopic;
                //    theoMasses[4, 0] = check5.Mass.Monoisotopic;
                //    theoMasses[5, 0] = check6.Mass.Monoisotopic;
                //    theoMasses[6, 0] = check7.Mass.Monoisotopic;
                //    theoMasses[7, 0] = check8.Mass.Monoisotopic;
                //    theoMasses[8, 0] = check9.Mass.Monoisotopic;
                //    theoMasses[9, 0] = check10.Mass.Monoisotopic;
                //    theoMasses[10, 0] = check11.Mass.Monoisotopic;
                //    theoMasses[11, 0] = check12.Mass.Monoisotopic;
                //}

                // 6plex NeuCode (chemical)
                //else if (Form1.NEUCODE_SIXPLEX_MTRAQ)
                //{
                //    check1 = new Peptide(peptide);
                //    check2 = new Peptide(peptide);
                //    check3 = new Peptide(peptide);
                //    check4 = new Peptide(peptide);
                //    check5 = new Peptide(peptide);
                //    check6 = new Peptide(peptide);

                //    if (tyrosineNHSPositions.Count > 0)
                //    {
                //        foreach (int position in tyrosineNHSPositions)
                //        {
                //            check1.SetModification(NamedChemicalFormula.GetModification("mTRAQ L"), position);
                //            check2.SetModification(NamedChemicalFormula.GetModification("mTRAQ L"), position);
                //            check3.SetModification(NamedChemicalFormula.GetModification("mTRAQ M"), position);
                //            check4.SetModification(NamedChemicalFormula.GetModification("mTRAQ M"), position);
                //            check5.SetModification(NamedChemicalFormula.GetModification("mTRAQ H"), position);
                //            check6.SetModification(NamedChemicalFormula.GetModification("mTRAQ H"), position);
                //        }
                //    }

                //    check1.SetModification(NamedChemicalFormula.GetModification("mTRAQ L"), ModificationSites.NPep);
                //    check2.SetModification(NamedChemicalFormula.GetModification("mTRAQ L"), ModificationSites.NPep);
                //    check3.SetModification(NamedChemicalFormula.GetModification("mTRAQ M"), ModificationSites.NPep);
                //    check4.SetModification(NamedChemicalFormula.GetModification("mTRAQ M"), ModificationSites.NPep);
                //    check5.SetModification(NamedChemicalFormula.GetModification("mTRAQ H"), ModificationSites.NPep);
                //    check6.SetModification(NamedChemicalFormula.GetModification("mTRAQ H"), ModificationSites.NPep);

                //    check1.SetModification(NamedChemicalFormula.GetModification("mTRAQ L Lys +8 13C6 15N2"), ModificationSites.K);
                //    check2.SetModification(NamedChemicalFormula.GetModification("mTRAQ L Lys +8 2H8"), ModificationSites.K);
                //    check3.SetModification(NamedChemicalFormula.GetModification("mTRAQ M Lys +8 13C6 15N2"), ModificationSites.K);
                //    check4.SetModification(NamedChemicalFormula.GetModification("mTRAQ M Lys +8 2H8"), ModificationSites.K);
                //    check5.SetModification(NamedChemicalFormula.GetModification("mTRAQ H Lys +8 13C6 15N2"), ModificationSites.K);
                //    check6.SetModification(NamedChemicalFormula.GetModification("mTRAQ H Lys +8 2H8"), ModificationSites.K);

                //    theoMasses[0, 0] = check1.Mass.Monoisotopic;
                //    theoMasses[1, 0] = check2.Mass.Monoisotopic;
                //    theoMasses[2, 0] = check3.Mass.Monoisotopic;
                //    theoMasses[3, 0] = check4.Mass.Monoisotopic;
                //    theoMasses[4, 0] = check5.Mass.Monoisotopic;
                //    theoMasses[5, 0] = check6.Mass.Monoisotopic;
                //}

                // 9plex NeuCode (chemical)
                //else if (Form1.NEUCODE_9PLEX_MTRAQ)
                //{
                //    check1 = new Peptide(peptide);
                //    check2 = new Peptide(peptide);
                //    check3 = new Peptide(peptide);
                //    check4 = new Peptide(peptide);
                //    check5 = new Peptide(peptide);
                //    check6 = new Peptide(peptide);
                //    check7 = new Peptide(peptide);
                //    check8 = new Peptide(peptide);
                //    check9 = new Peptide(peptide);

                //    if (tyrosineNHSPositions.Count > 0)
                //    {
                //        foreach (int position in tyrosineNHSPositions)
                //        {
                //            check1.SetModification(NamedChemicalFormula.GetModification("mTRAQ L"), position);
                //            check2.SetModification(NamedChemicalFormula.GetModification("mTRAQ L"), position);
                //            check3.SetModification(NamedChemicalFormula.GetModification("mTRAQ L"), position);
                //            check4.SetModification(NamedChemicalFormula.GetModification("mTRAQ M"), position);
                //            check5.SetModification(NamedChemicalFormula.GetModification("mTRAQ M"), position);
                //            check6.SetModification(NamedChemicalFormula.GetModification("mTRAQ M"), position);
                //            check7.SetModification(NamedChemicalFormula.GetModification("mTRAQ H"), position);
                //            check8.SetModification(NamedChemicalFormula.GetModification("mTRAQ H"), position);
                //            check9.SetModification(NamedChemicalFormula.GetModification("mTRAQ H"), position);
                //        }
                //    }

                //    check1.SetModification(NamedChemicalFormula.GetModification("mTRAQ L"), ModificationSites.NPep);
                //    check2.SetModification(NamedChemicalFormula.GetModification("mTRAQ L"), ModificationSites.NPep);
                //    check3.SetModification(NamedChemicalFormula.GetModification("mTRAQ L"), ModificationSites.NPep);
                //    check4.SetModification(NamedChemicalFormula.GetModification("mTRAQ M"), ModificationSites.NPep);
                //    check5.SetModification(NamedChemicalFormula.GetModification("mTRAQ M"), ModificationSites.NPep);
                //    check6.SetModification(NamedChemicalFormula.GetModification("mTRAQ M"), ModificationSites.NPep);
                //    check7.SetModification(NamedChemicalFormula.GetModification("mTRAQ H"), ModificationSites.NPep);
                //    check8.SetModification(NamedChemicalFormula.GetModification("mTRAQ H"), ModificationSites.NPep);
                //    check9.SetModification(NamedChemicalFormula.GetModification("mTRAQ H"), ModificationSites.NPep);

                //    check1.SetModification(NamedChemicalFormula.GetModification("mTRAQ L Lys +8 13C6 15N2"), ModificationSites.K);
                //    check2.SetModification(NamedChemicalFormula.GetModification("mTRAQ L Lys +8 2H6 15N2"), ModificationSites.K);
                //    check3.SetModification(NamedChemicalFormula.GetModification("mTRAQ L Lys +8 2H8"), ModificationSites.K);
                //    check4.SetModification(NamedChemicalFormula.GetModification("mTRAQ M Lys +8 13C6 15N2"), ModificationSites.K);
                //    check5.SetModification(NamedChemicalFormula.GetModification("mTRAQ M Lys +8 2H6 15N2"), ModificationSites.K);
                //    check6.SetModification(NamedChemicalFormula.GetModification("mTRAQ M Lys +8 2H8"), ModificationSites.K);
                //    check7.SetModification(NamedChemicalFormula.GetModification("mTRAQ H Lys +8 13C6 15N2"), ModificationSites.K);
                //    check8.SetModification(NamedChemicalFormula.GetModification("mTRAQ H Lys +8 2H6 15N2"), ModificationSites.K);
                //    check9.SetModification(NamedChemicalFormula.GetModification("mTRAQ H Lys +8 2H8"), ModificationSites.K);

                //    theoMasses[0, 0] = check1.Mass.Monoisotopic;
                //    theoMasses[1, 0] = check2.Mass.Monoisotopic;
                //    theoMasses[2, 0] = check3.Mass.Monoisotopic;
                //    theoMasses[3, 0] = check4.Mass.Monoisotopic;
                //    theoMasses[4, 0] = check5.Mass.Monoisotopic;
                //    theoMasses[5, 0] = check6.Mass.Monoisotopic;
                //    theoMasses[6, 0] = check7.Mass.Monoisotopic;
                //    theoMasses[7, 0] = check8.Mass.Monoisotopic;
                //    theoMasses[8, 0] = check9.Mass.Monoisotopic;
                //}

                // 12plex NeuCode (chemical)
                //else if (Form1.NEUCODE_12PLEX_MTRAQ)
                //{
                //    check1 = new Peptide(peptide);
                //    check2 = new Peptide(peptide);
                //    check3 = new Peptide(peptide);
                //    check4 = new Peptide(peptide);
                //    check5 = new Peptide(peptide);
                //    check6 = new Peptide(peptide);
                //    check7 = new Peptide(peptide);
                //    check8 = new Peptide(peptide);
                //    check9 = new Peptide(peptide);
                //    check10 = new Peptide(peptide);
                //    check11 = new Peptide(peptide);
                //    check12 = new Peptide(peptide);

                //    if (tyrosineNHSPositions.Count > 0)
                //    {
                //        foreach (int position in tyrosineNHSPositions)
                //        {
                //            check1.SetModification(NamedChemicalFormula.GetModification("mTRAQ L"), position);
                //            check2.SetModification(NamedChemicalFormula.GetModification("mTRAQ L"), position);
                //            check3.SetModification(NamedChemicalFormula.GetModification("mTRAQ L"), position);
                //            check4.SetModification(NamedChemicalFormula.GetModification("mTRAQ L"), position);
                //            check5.SetModification(NamedChemicalFormula.GetModification("mTRAQ M"), position);
                //            check6.SetModification(NamedChemicalFormula.GetModification("mTRAQ M"), position);
                //            check7.SetModification(NamedChemicalFormula.GetModification("mTRAQ M"), position);
                //            check8.SetModification(NamedChemicalFormula.GetModification("mTRAQ M"), position);
                //            check9.SetModification(NamedChemicalFormula.GetModification("mTRAQ H"), position);
                //            check10.SetModification(NamedChemicalFormula.GetModification("mTRAQ H"), position);
                //            check11.SetModification(NamedChemicalFormula.GetModification("mTRAQ H"), position);
                //            check12.SetModification(NamedChemicalFormula.GetModification("mTRAQ H"), position);
                //        }
                //    }

                //    check1.SetModification(NamedChemicalFormula.GetModification("mTRAQ L"), ModificationSites.NPep);
                //    check2.SetModification(NamedChemicalFormula.GetModification("mTRAQ L"), ModificationSites.NPep);
                //    check3.SetModification(NamedChemicalFormula.GetModification("mTRAQ L"), ModificationSites.NPep);
                //    check4.SetModification(NamedChemicalFormula.GetModification("mTRAQ L"), ModificationSites.NPep);
                //    check5.SetModification(NamedChemicalFormula.GetModification("mTRAQ M"), ModificationSites.NPep);
                //    check6.SetModification(NamedChemicalFormula.GetModification("mTRAQ M"), ModificationSites.NPep);
                //    check7.SetModification(NamedChemicalFormula.GetModification("mTRAQ M"), ModificationSites.NPep);
                //    check8.SetModification(NamedChemicalFormula.GetModification("mTRAQ M"), ModificationSites.NPep);
                //    check9.SetModification(NamedChemicalFormula.GetModification("mTRAQ H"), ModificationSites.NPep);
                //    check10.SetModification(NamedChemicalFormula.GetModification("mTRAQ H"), ModificationSites.NPep);
                //    check11.SetModification(NamedChemicalFormula.GetModification("mTRAQ H"), ModificationSites.NPep);
                //    check12.SetModification(NamedChemicalFormula.GetModification("mTRAQ H"), ModificationSites.NPep);

                //    check1.SetModification(NamedChemicalFormula.GetModification("mTRAQ L Lys +8 13C6 15N2"), ModificationSites.K);
                //    check2.SetModification(NamedChemicalFormula.GetModification("mTRAQ L Lys +8 13C5 2H2 15N1"), ModificationSites.K);
                //    check3.SetModification(NamedChemicalFormula.GetModification("mTRAQ L Lys +8 13C4 2H4"), ModificationSites.K);
                //    check4.SetModification(NamedChemicalFormula.GetModification("mTRAQ L Lys +8 2H8"), ModificationSites.K);
                //    check5.SetModification(NamedChemicalFormula.GetModification("mTRAQ M Lys +8 13C6 15N2"), ModificationSites.K);
                //    check6.SetModification(NamedChemicalFormula.GetModification("mTRAQ M Lys +8 13C5 2H2 15N1"), ModificationSites.K);
                //    check7.SetModification(NamedChemicalFormula.GetModification("mTRAQ M Lys +8 13C4 2H4"), ModificationSites.K);
                //    check8.SetModification(NamedChemicalFormula.GetModification("mTRAQ M Lys +8 2H8"), ModificationSites.K);
                //    check9.SetModification(NamedChemicalFormula.GetModification("mTRAQ H Lys +8 13C6 15N2"), ModificationSites.K);
                //    check10.SetModification(NamedChemicalFormula.GetModification("mTRAQ H Lys +8 13C5 2H2 15N1"), ModificationSites.K);
                //    check11.SetModification(NamedChemicalFormula.GetModification("mTRAQ H Lys +8 13C4 2H4"), ModificationSites.K);
                //    check12.SetModification(NamedChemicalFormula.GetModification("mTRAQ H Lys +8 2H8"), ModificationSites.K);

                //    theoMasses[0, 0] = check1.Mass.Monoisotopic;
                //    theoMasses[1, 0] = check2.Mass.Monoisotopic;
                //    theoMasses[2, 0] = check3.Mass.Monoisotopic;
                //    theoMasses[3, 0] = check4.Mass.Monoisotopic;
                //    theoMasses[4, 0] = check5.Mass.Monoisotopic;
                //    theoMasses[5, 0] = check6.Mass.Monoisotopic;
                //    theoMasses[6, 0] = check7.Mass.Monoisotopic;
                //    theoMasses[7, 0] = check8.Mass.Monoisotopic;
                //    theoMasses[8, 0] = check9.Mass.Monoisotopic;
                //    theoMasses[9, 0] = check10.Mass.Monoisotopic;
                //    theoMasses[10, 0] = check11.Mass.Monoisotopic;
                //    theoMasses[11, 0] = check12.Mass.Monoisotopic;
                //}

                // 18plex NeuCode (chemical)
                //else if (Form1.NEUCODE_18PLEX_MTRAQ)
                //{
                //    check1 = new Peptide(peptide);
                //    check2 = new Peptide(peptide);
                //    check3 = new Peptide(peptide);
                //    check4 = new Peptide(peptide);
                //    check5 = new Peptide(peptide);
                //    check6 = new Peptide(peptide);
                //    check7 = new Peptide(peptide);
                //    check8 = new Peptide(peptide);
                //    check9 = new Peptide(peptide);
                //    check10 = new Peptide(peptide);
                //    check11 = new Peptide(peptide);
                //    check12 = new Peptide(peptide);
                //    check13 = new Peptide(peptide);
                //    check14 = new Peptide(peptide);
                //    check15 = new Peptide(peptide);
                //    check16 = new Peptide(peptide);
                //    check17 = new Peptide(peptide);
                //    check18 = new Peptide(peptide);

                //    if (tyrosineNHSPositions.Count > 0)
                //    {
                //        foreach (int position in tyrosineNHSPositions)
                //        {
                //            check1.SetModification(NamedChemicalFormula.GetModification("mTRAQ L"), position);
                //            check2.SetModification(NamedChemicalFormula.GetModification("mTRAQ L"), position);
                //            check3.SetModification(NamedChemicalFormula.GetModification("mTRAQ L"), position);
                //            check4.SetModification(NamedChemicalFormula.GetModification("mTRAQ L"), position);
                //            check5.SetModification(NamedChemicalFormula.GetModification("mTRAQ L"), position);
                //            check6.SetModification(NamedChemicalFormula.GetModification("mTRAQ L"), position);
                //            check7.SetModification(NamedChemicalFormula.GetModification("mTRAQ M"), position);
                //            check8.SetModification(NamedChemicalFormula.GetModification("mTRAQ M"), position);
                //            check9.SetModification(NamedChemicalFormula.GetModification("mTRAQ M"), position);
                //            check10.SetModification(NamedChemicalFormula.GetModification("mTRAQ M"), position);
                //            check11.SetModification(NamedChemicalFormula.GetModification("mTRAQ M"), position);
                //            check12.SetModification(NamedChemicalFormula.GetModification("mTRAQ M"), position);
                //            check13.SetModification(NamedChemicalFormula.GetModification("mTRAQ H"), position);
                //            check14.SetModification(NamedChemicalFormula.GetModification("mTRAQ H"), position);
                //            check15.SetModification(NamedChemicalFormula.GetModification("mTRAQ H"), position);
                //            check16.SetModification(NamedChemicalFormula.GetModification("mTRAQ H"), position);
                //            check17.SetModification(NamedChemicalFormula.GetModification("mTRAQ H"), position);
                //            check18.SetModification(NamedChemicalFormula.GetModification("mTRAQ H"), position);
                //        }
                //    }

                //    check1.SetModification(NamedChemicalFormula.GetModification("mTRAQ L"), ModificationSites.NPep);
                //    check2.SetModification(NamedChemicalFormula.GetModification("mTRAQ L"), ModificationSites.NPep);
                //    check3.SetModification(NamedChemicalFormula.GetModification("mTRAQ L"), ModificationSites.NPep);
                //    check4.SetModification(NamedChemicalFormula.GetModification("mTRAQ L"), ModificationSites.NPep);
                //    check5.SetModification(NamedChemicalFormula.GetModification("mTRAQ L"), ModificationSites.NPep);
                //    check6.SetModification(NamedChemicalFormula.GetModification("mTRAQ L"), ModificationSites.NPep);
                //    check7.SetModification(NamedChemicalFormula.GetModification("mTRAQ M"), ModificationSites.NPep);
                //    check8.SetModification(NamedChemicalFormula.GetModification("mTRAQ M"), ModificationSites.NPep);
                //    check9.SetModification(NamedChemicalFormula.GetModification("mTRAQ M"), ModificationSites.NPep);
                //    check10.SetModification(NamedChemicalFormula.GetModification("mTRAQ M"), ModificationSites.NPep);
                //    check11.SetModification(NamedChemicalFormula.GetModification("mTRAQ M"), ModificationSites.NPep);
                //    check12.SetModification(NamedChemicalFormula.GetModification("mTRAQ M"), ModificationSites.NPep);
                //    check13.SetModification(NamedChemicalFormula.GetModification("mTRAQ H"), ModificationSites.NPep);
                //    check14.SetModification(NamedChemicalFormula.GetModification("mTRAQ H"), ModificationSites.NPep);
                //    check15.SetModification(NamedChemicalFormula.GetModification("mTRAQ H"), ModificationSites.NPep);
                //    check16.SetModification(NamedChemicalFormula.GetModification("mTRAQ H"), ModificationSites.NPep);
                //    check17.SetModification(NamedChemicalFormula.GetModification("mTRAQ H"), ModificationSites.NPep);
                //    check18.SetModification(NamedChemicalFormula.GetModification("mTRAQ H"), ModificationSites.NPep);

                //    check1.SetModification(NamedChemicalFormula.GetModification("mTRAQ L Lys +8 13C6 15N2"), ModificationSites.K);
                //    check2.SetModification(NamedChemicalFormula.GetModification("mTRAQ L Lys +8 13C4 2H2 15N2"), ModificationSites.K);
                //    check3.SetModification(NamedChemicalFormula.GetModification("mTRAQ L Lys +8 13C5 2H2 15N1"), ModificationSites.K);
                //    check4.SetModification(NamedChemicalFormula.GetModification("mTRAQ L Lys +8 2H6 15N2"), ModificationSites.K);
                //    check5.SetModification(NamedChemicalFormula.GetModification("mTRAQ L Lys +8 13C4 2H4"), ModificationSites.K);
                //    check6.SetModification(NamedChemicalFormula.GetModification("mTRAQ L Lys +8 2H8"), ModificationSites.K);
                //    check7.SetModification(NamedChemicalFormula.GetModification("mTRAQ M Lys +8 13C6 15N2"), ModificationSites.K);
                //    check8.SetModification(NamedChemicalFormula.GetModification("mTRAQ M Lys +8 13C4 2H2 15N2"), ModificationSites.K);
                //    check9.SetModification(NamedChemicalFormula.GetModification("mTRAQ M Lys +8 13C5 2H2 15N1"), ModificationSites.K);
                //    check10.SetModification(NamedChemicalFormula.GetModification("mTRAQ M Lys +8 2H6 15N2"), ModificationSites.K);
                //    check11.SetModification(NamedChemicalFormula.GetModification("mTRAQ M Lys +8 13C4 2H4"), ModificationSites.K);
                //    check12.SetModification(NamedChemicalFormula.GetModification("mTRAQ M Lys +8 2H8"), ModificationSites.K);
                //    check13.SetModification(NamedChemicalFormula.GetModification("mTRAQ H Lys +8 13C6 15N2"), ModificationSites.K);
                //    check14.SetModification(NamedChemicalFormula.GetModification("mTRAQ H Lys +8 13C4 2H2 15N2"), ModificationSites.K);
                //    check15.SetModification(NamedChemicalFormula.GetModification("mTRAQ H Lys +8 13C5 2H2 15N1"), ModificationSites.K);
                //    check16.SetModification(NamedChemicalFormula.GetModification("mTRAQ H Lys +8 2H6 15N2"), ModificationSites.K);
                //    check17.SetModification(NamedChemicalFormula.GetModification("mTRAQ H Lys +8 13C4 2H4"), ModificationSites.K);
                //    check18.SetModification(NamedChemicalFormula.GetModification("mTRAQ H Lys +8 2H8"), ModificationSites.K);

                //    theoMasses[0, 0] = check1.Mass.Monoisotopic;
                //    theoMasses[1, 0] = check2.Mass.Monoisotopic;
                //    theoMasses[2, 0] = check3.Mass.Monoisotopic;
                //    theoMasses[3, 0] = check4.Mass.Monoisotopic;
                //    theoMasses[4, 0] = check5.Mass.Monoisotopic;
                //    theoMasses[5, 0] = check6.Mass.Monoisotopic;
                //    theoMasses[6, 0] = check7.Mass.Monoisotopic;
                //    theoMasses[7, 0] = check8.Mass.Monoisotopic;
                //    theoMasses[8, 0] = check9.Mass.Monoisotopic;
                //    theoMasses[9, 0] = check10.Mass.Monoisotopic;
                //    theoMasses[10, 0] = check11.Mass.Monoisotopic;
                //    theoMasses[11, 0] = check12.Mass.Monoisotopic;
                //    theoMasses[12, 0] = check13.Mass.Monoisotopic;
                //    theoMasses[13, 0] = check14.Mass.Monoisotopic;
                //    theoMasses[14, 0] = check15.Mass.Monoisotopic;
                //    theoMasses[15, 0] = check16.Mass.Monoisotopic;
                //    theoMasses[16, 0] = check17.Mass.Monoisotopic;
                //    theoMasses[17, 0] = check18.Mass.Monoisotopic;
                //}

                // 6plex NeuCode (metabolic)
                //else if (Form1.NEUCODE_SIXPLEX_ARG || Form1.NEUCODE_SIXPLEX_LEU)
                //{
                //    check1 = new Peptide(peptide);
                //    check2 = new Peptide(peptide);
                //    check3 = new Peptide(peptide);
                //    check4 = new Peptide(peptide);
                //    check5 = new Peptide(peptide);
                //    check6 = new Peptide(peptide);

                //    check1.SetModification(NamedChemicalFormula.GetModification("Lys +8 13C6 15N2"), ModificationSites.K);
                //    check2.SetModification(NamedChemicalFormula.GetModification("Lys +8 2H8"), ModificationSites.K);

                //    if ((Form1.NEUCODE_SIXPLEX_ARG && countResidues('R', peptide.Sequence) < 1) || Form1.NEUCODE_SIXPLEX_LEU && countResidues('L', peptide.Sequence) < 1)
                //    {
                //        theoMasses[0, 0] = check1.Mass.Monoisotopic;
                //        theoMasses[1, 0] = check2.Mass.Monoisotopic;
                //    }
                //    else
                //    {
                //        check3.SetModification(NamedChemicalFormula.GetModification("Lys +8 13C6 15N2"), ModificationSites.K);
                //        check4.SetModification(NamedChemicalFormula.GetModification("Lys +8 2H8"), ModificationSites.K);
                //        check5.SetModification(NamedChemicalFormula.GetModification("Lys +8 13C6 15N2"), ModificationSites.K);
                //        check6.SetModification(NamedChemicalFormula.GetModification("Lys +8 2H8"), ModificationSites.K);

                //        if (Form1.NEUCODE_SIXPLEX_ARG)
                //        {
                //            check3.SetModification(NamedChemicalFormula.GetModification("Arg +6 13C6"), ModificationSites.R);
                //            check4.SetModification(NamedChemicalFormula.GetModification("Arg +6 13C6"), ModificationSites.R);
                //            check5.SetModification(NamedChemicalFormula.GetModification("Arg +10 13C6 15N4"), ModificationSites.R);
                //            check6.SetModification(NamedChemicalFormula.GetModification("Arg +10 13C6 15N4"), ModificationSites.R);
                //        }
                //        else
                //        {
                //            check3.SetModification(NamedChemicalFormula.GetModification("Leu +6 13C6 15N1"), ModificationSites.L);
                //            check4.SetModification(NamedChemicalFormula.GetModification("Leu +6 13C6 15N1"), ModificationSites.L);
                //            check5.SetModification(NamedChemicalFormula.GetModification("Leu +10 2H10"), ModificationSites.L);
                //            check6.SetModification(NamedChemicalFormula.GetModification("Leu +10 2H10"), ModificationSites.L);
                //        }

                //        theoMasses[0, 0] = check1.Mass.Monoisotopic;
                //        theoMasses[1, 0] = check2.Mass.Monoisotopic;
                //        theoMasses[2, 0] = check3.Mass.Monoisotopic;
                //        theoMasses[3, 0] = check4.Mass.Monoisotopic;
                //        theoMasses[4, 0] = check5.Mass.Monoisotopic;
                //        theoMasses[5, 0] = check6.Mass.Monoisotopic;
                //    }
                //}

                // Duplex SILAC
                //else if (Form1.SILAC_DUPLEX_LYSC || Form1.SILAC_DUPLEX_LYSCN || Form1.SILAC_DUPLEX_LYSH || Form1.SILAC_DUPLEX_LEUCN || Form1.SILAC_DUPLEX_LEUH)
                //{
                //    check1 = new Peptide(peptide);
                //    check2 = new Peptide(peptide);

                //    if (Form1.SILAC_DUPLEX_LYSC)
                //    {
                //        check2.SetModification(NamedChemicalFormula.GetModification("Lys +6 13C6"), ModificationSites.K);
                //    }
                //    else if (Form1.SILAC_DUPLEX_LYSCN)
                //    {
                //        check2.SetModification(NamedChemicalFormula.GetModification("Lys +8 13C6 15N2"), ModificationSites.K);
                //    }
                //    else if (Form1.SILAC_DUPLEX_LYSH)
                //    {
                //        check2.SetModification(NamedChemicalFormula.GetModification("Lys +8 2H8"), ModificationSites.K);
                //    }
                //    else if (Form1.SILAC_DUPLEX_LEUCN)
                //    {
                //        check2.SetModification(NamedChemicalFormula.GetModification("Leu +7 13C6 15N1"), ModificationSites.L);
                //    }
                //    else
                //    {
                //        check2.SetModification(NamedChemicalFormula.GetModification("Leu +7 2H7"), ModificationSites.L);
                //    }

                //    theoMasses[0, 0] = check1.Mass.Monoisotopic;
                //    theoMasses[1, 0] = check2.Mass.Monoisotopic;
                //}
                //else if (Form1.ICAT)
                //{
                //    check1 = new Peptide(peptide);
                //    check2 = new Peptide(peptide);

                //    check1.SetModification(NamedChemicalFormula.GetModification("iCAT L"), ModificationSites.C);
                //    check2.SetModification(NamedChemicalFormula.GetModification("iCAT H"), ModificationSites.C);

                //    theoMasses[0, 0] = check1.Mass.Monoisotopic;
                //    theoMasses[1, 0] = check2.Mass.Monoisotopic;
                //}
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
                    //int peakCharge = 0;
                    //// Ignore any peak with a known charge state that is incorrect
                    //if (largest.Charge != charge)
                    //{
                    //    MSDataScan previous = rawFile[currentScan.SpectrumNumber - 1];
                    //    peakCharge = findLowerResolutionCharge(largest, previous);
                    //}
                    
                    if (largest.Charge > 0 && largest.Charge != charge) 
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

        public void maximizeResolvability()
        {
            foreach (KeyValuePair<MSDataFile, List<PeptideSpectralMatch>> psmSet in PSMs)
            {
                List<PeptideSpectralMatch> psms = psmSet.Value;
                PeptideSpectralMatch bestResolvablePSM;
                double bestTheoResolvability;
                double bestExpResolvability;
                double bestResolvability;
                PeptideSpectralMatch currentPSM;
                double currentTheoResolvability;
                double currentExpResolvability;
                double currentResolvability;

                if (psms.Count > 1 && numIsotopologueLabels > 0)
                {
                    double coefficient = (Math.Sqrt(2 * Math.Log(100.0 / QUANTSEPARATION))) / (Math.Sqrt(2 * Math.Log(2)));
                    double spacing = spacingMassRange[1,0].Mean;
                    bestResolvablePSM = psms[0];
                    bestTheoResolvability = coefficient * ((Mass.MzFromMass(theoMasses[0, 0], bestResolvablePSM.Charge) / (QUANTRESOLUTION * Math.Sqrt(400 / Mass.MzFromMass(theoMasses[0, 0], bestResolvablePSM.Charge)))));
                    bestExpResolvability = spacing / (double)bestResolvablePSM.Charge;
                    bestResolvability = bestExpResolvability / bestTheoResolvability;

                    for (int i = 1; i < psms.Count; i++)
                    {
                        currentPSM = psms[i];
                        currentTheoResolvability = coefficient * ((Mass.MzFromMass(theoMasses[0, 0], currentPSM.Charge) / (QUANTRESOLUTION * Math.Sqrt(400 / Mass.MzFromMass(theoMasses[0, 0], currentPSM.Charge)))));
                        currentExpResolvability = spacing / (double)currentPSM.Charge;
                        currentResolvability = currentExpResolvability / currentTheoResolvability;
                        if (currentResolvability > bestResolvability)
                        {
                            bestResolvablePSM = currentPSM;
                            bestResolvability = currentResolvability;
                        }
                    }

                    bestPSMs[psmSet.Key] = bestResolvablePSM;
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
            double resolutionMinimum = 100000.0;

            if (Form1.FUSION)
            {
                resolutionMinimum = 0.0;
            }
            else if (numIsotopologues < 2)
            {
                resolutionMinimum = 15000.0;
            }

            fullScanList = new List<MSDataScan>();

            // For peptides eluting across more than 1000 scans, create retention time window around best PSM
            if (last - first < 1000)
            {
                rTRange = new Range<double>(firstScanTime - rTWindow, lastScanTime + (2 * rTWindow));
            }
            else
            {
                rTRange = new Range<double>(bestScanTime - (2 * rTWindow), bestScanTime + (2 * rTWindow));
            }

            // Check to see if retention time window extends past the length of the run
            if (rTRange.Minimum < firstRunTime)
            {
                Range<double> lowExceeded = new Range<double>(firstRunTime, rTRange.Maximum);
                rTRange = lowExceeded; //Minimum time surpassed
            }

            if (rTRange.Maximum > lastRunTime)
            {
                Range<double> highExceeded = new Range<double>(rTRange.Minimum, lastRunTime);
                rTRange = highExceeded; //Maximum time surpassed
            }

            int firstScanNumber = rawFile.GetSpectrumNumber(rTRange.Minimum);
            int lastScanNumber = rawFile.GetSpectrumNumber(rTRange.Maximum);

            // For NeuCode, use high resolution MS1 scans
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
            //else if (Form1.AGCBINNING)
            //{
            //    foreach (MSDataScan scan in rawFile.Where(scan => scan.MsnOrder == 2 && ((MsnDataScan)scan).DissociationType == DissociationType.HCD && scan.RetentionTime <= rTRange.Maximum && scan.RetentionTime >= rTRange.Minimum))
            //    {
            //        fullScanList.Add(scan);
            //    }
            //}
            //else if (Form1.NEUCODE || Form1.SILAC_DUPLEX_LEUCN || Form1.SILAC_DUPLEX_LEUH || Form1.SILAC_DUPLEX_LYSCN || Form1.SILAC_DUPLEX_LYSH)
            //{
            //    foreach (MSDataScan scan in rawFile.Where(scan => scan.MsnOrder == 1 && scan.Resolution >= resolutionMinimum && scan.RetentionTime <= rTRange.Maximum && scan.RetentionTime >= rTRange.Minimum))
            //    {
            //        fullScanList.Add(scan);
            //    }
            //}
            //Otherwise, include every MS1 scan for quantification
            else
            {
                foreach (MSDataScan scan in rawFile.Where(scan => scan.MsnOrder == 1 && scan.Resolution >= resolutionMinimum && scan.RetentionTime <= rTRange.Maximum && scan.RetentionTime >= rTRange.Minimum))
                {
                    fullScanList.Add(scan);
                }
            }
        }

        public List<ILabeledPeak> findPeakPatterns(MassRange Range, int Charge, MSDataScan CurrentScan, List<Spacing> Spacings = null, int Isotope = 0)
        {
            List<MZPeak> patternPeaksFound = new List<MZPeak>();
            List<ILabeledPeak> convertedPeaks = null;
            if (CurrentScan.MassSpectrum.TryGetPeaks(Range.Minimum, Range.Maximum, out patternPeaksFound))
            {
                convertedPeaks = new List<ILabeledPeak>();

                foreach (IPeak peakToConvert in patternPeaksFound)
                {
                    ILabeledPeak newPeak = (ILabeledPeak)peakToConvert;
                    if (newPeak.GetSignalToNoise() < SIGNALTONOISE)
                    {
                        newPeak = null;
                    }
                    else if (newPeak.Charge != 0 && newPeak.Charge != Charge)
                    {
                        newPeak = null;
                    }
                    else
                    {
                        convertedPeaks.Add(newPeak);
                    }
                }

                if (convertedPeaks.Count > 0 && Form1.OUTPUTSPACINGS)
                {
                    List<ILabeledPeak> intensitySortedPeaks = convertedPeaks.OrderByDescending(peak => peak.GetSignalToNoise()).Take(2).ToList();
                    List<ILabeledPeak> mZSortedPeaks = intensitySortedPeaks.OrderBy(peak => peak.X).ToList();

                    if (mZSortedPeaks.Count == 1)
                    {
                        Spacings.Add(new Spacing(sequence, CurrentScan.RetentionTime, Charge, numIsotopologueLabels, Isotope, this, mZSortedPeaks[0]));
                    }
                    else if (mZSortedPeaks.Count == 2)
                    {
                        Spacings.Add(new Spacing(sequence, CurrentScan.RetentionTime, Charge, numIsotopologueLabels, Isotope, this, mZSortedPeaks[0], mZSortedPeaks[1]));
                    }                    
                }

                convertedPeaks = cleanPeaks(convertedPeaks, Charge);
            }
            return convertedPeaks;
        }

        public MassRange checkMappedPeaks(ILabeledPeak[] MappedPeaks, double ToleranceWidth, MassTolerance Tolerance)
        {
            MassRange newSearchWindow = null;
            double halfWidth = ToleranceWidth / 2.0;
            List<int> peakPositions = new List<int>();

            for (int i = 0; i < MappedPeaks.Length; i++)
            {
                if (MappedPeaks[i] != null && MappedPeaks[i].GetSignalToNoise() >= Form1.MINIMUMSN)
                {
                    peakPositions.Add(i);
                }
            }

            int peakCount = peakPositions.Count;

            if (peakCount >= peaksNeeded)
            {
                double newMinMZ = 0;
                double newMaxMZ = 0;
                if (numIsotopologues == 3)
                {
                    if (peakCount == 1)
                    {
                        if (peakPositions[0] == 0)
                        {
                            newMaxMZ = new MassRange(MappedPeaks[0].X, Tolerance).Maximum;
                            newMinMZ = newMaxMZ - ToleranceWidth;
                        }
                        else if (peakPositions[0] == 2)
                        {
                            newMinMZ = new MassRange(MappedPeaks[2].X, Tolerance).Minimum;
                            newMaxMZ = newMinMZ + ToleranceWidth;
                        }
                    }
                    else if (peakCount == 2)
                    {
                        if (peakPositions[0] == 0 && peakPositions[1] == 1)
                        {
                            newMaxMZ = new MassRange(MappedPeaks[1].X, Tolerance).Maximum;
                            newMinMZ = newMaxMZ - ToleranceWidth;
                        }
                        else if (peakPositions[0] == 1 && peakPositions[1] == 2)
                        {
                            newMinMZ = new MassRange(MappedPeaks[1].X, Tolerance).Minimum;
                            newMaxMZ = newMinMZ + ToleranceWidth;
                        }
                    }
                }
                if (numIsotopologues == 4)
                {
                    if (peakCount == 2)
                    {
                        if (peakPositions[0] == 0 && peakPositions[1] == 1)
                        {
                            newMaxMZ = new MassRange(MappedPeaks[1].X, Tolerance).Maximum;
                            newMinMZ = newMaxMZ - ToleranceWidth;
                        }
                        else if(peakPositions[0] == 0 && peakPositions[1] == 2)
                        {
                            newMaxMZ = new MassRange(MappedPeaks[2].X, Tolerance).Maximum;
                            newMinMZ = newMaxMZ - ToleranceWidth;
                        }
                        else if (peakPositions[0] == 2 && peakPositions[1] == 3)
                        {
                            newMinMZ = new MassRange(MappedPeaks[2].X, Tolerance).Minimum;
                            newMaxMZ = newMinMZ + ToleranceWidth;
                        }
                        else if (peakPositions[0] == 1 && peakPositions[1] == 3)
                        {
                            newMinMZ = new MassRange(MappedPeaks[1].X, Tolerance).Minimum;
                            newMaxMZ = newMinMZ + ToleranceWidth;
                        }
                    }
                    else if (peakCount == 3)
                    {
                        if (peakPositions[0] == 0 && peakPositions[1] == 1 && peakPositions[2] == 2)
                        {
                            newMaxMZ = new MassRange(MappedPeaks[2].X, Tolerance).Maximum;
                            newMinMZ = newMaxMZ - ToleranceWidth;
                        }
                        else if (peakPositions[0] == 1 && peakPositions[1] == 2 && peakPositions[2] == 3)
                        {
                            newMinMZ = new MassRange(MappedPeaks[1].X, Tolerance).Minimum;
                            newMaxMZ = newMinMZ + ToleranceWidth;
                        }
                    }
                }
                newSearchWindow = new MassRange(newMinMZ, newMaxMZ);                
            }

            return newSearchWindow;
        }
        
        /* Finds light and heavy peaks in the current scan, using a specified PPM window centered around the peptide's adjusted masses
         * Adds a pair to the peptide's list of all pairs for that raw file if enough (i.e., equal to or greater than the number of channels) peaks are found
         */
        public void findPeaks(MSDataScan current, MSDataFile rawFile, List<Spacing> spacings = null)
        {
            int charge = bestPSMs[rawFile].Charge;
            ILabeledPeak[,] peaksFound = new ILabeledPeak[numChannels, numIsotopes];
            Pair pair;
            double injectionTime = rawFile.GetInjectionTime(current.SpectrumNumber);
            int peaksCount = 0;
            int isotopeCount = 0;
            int clusterCount = 0;
            bool stop = false;

            if (numIsotopologues > 1)
            {
                int channelIndex;

                clusterCount = 0;
                for (int c = 0; c < numClusters; c++) // Start cluster loop
                {
                    if (GetTheoreticalResolvability(c))
                    {
                        isotopeCount = 0;
                        for (int j = 0; j < numIsotopes; j++) // Start isotope loop
                        {
                            peaksCount = 0;
                            channelIndex = c * numIsotopologues;

                            // Start pattern search
                            MassTolerance tolerance = new MassTolerance(MassToleranceType.PPM, TOLERANCE.Value);

                            MassRange minRange = new MassRange(adjustedTheoMasses[channelIndex,j], tolerance);
                            MassRange maxRange = new MassRange(adjustedTheoMasses[channelIndex + (numIsotopologues - 1), j], tolerance);

                            double minMZ = Mass.MzFromMass(minRange.Minimum, charge);
                            double maxMZ = Mass.MzFromMass(maxRange.Maximum, charge);

                            MassRange totalRange = new MassRange(minMZ, maxMZ);

                            List<ILabeledPeak> cleanedPeaks = findPeakPatterns(totalRange, charge, current, spacings, j);
                            List<ILabeledPeak> topSignalToNoisePeaks = null;
                            List<ILabeledPeak> orderedByMzPeaks = null;
                            if (cleanedPeaks != null && cleanedPeaks.Count > 0)
                            {
                                topSignalToNoisePeaks = cleanedPeaks.OrderByDescending(peakSort => peakSort.GetSignalToNoise()).Take(numIsotopologues).ToList();
                                // Peaks found for a complete set of isotopologues                             
                                if (topSignalToNoisePeaks.Count == numIsotopologues)
                                {
                                    orderedByMzPeaks = topSignalToNoisePeaks.OrderBy(peakSort => peakSort.X).ToList();
                                    for (int i = 0; i < numIsotopologues; i++)
                                    {
                                        peaksFound[channelIndex + i, j] = orderedByMzPeaks[i];
                                    }
                                }
                                else
                                {
                                    if (!Form1.NOISEBANDCAP)
                                    {
                                        if (topSignalToNoisePeaks.Count >= numIsotopologues / 2)
                                        {
                                            ILabeledPeak[] mappedPeaks = mapPeaks(topSignalToNoisePeaks, channelIndex, channelIndex + (numIsotopologues - 1), j, charge);
                                            MassRange newSearchWindow = checkMappedPeaks(mappedPeaks, totalRange.Width, tolerance);

                                            if (newSearchWindow != null)
                                            {
                                                List<ILabeledPeak> secondTryPeaks = findPeakPatterns(newSearchWindow, charge, current);
                                                List<ILabeledPeak> secondTryTopSignalToNoisePeaks = null;
                                                List<ILabeledPeak> secondTryOrderedByMzPeaks = null;
                                                if (secondTryPeaks != null && secondTryPeaks.Count == numIsotopologues)
                                                {
                                                    secondTryTopSignalToNoisePeaks = secondTryPeaks.OrderByDescending(peakSort => peakSort.GetSignalToNoise()).Take(numIsotopologues).ToList();
                                                    secondTryOrderedByMzPeaks = secondTryTopSignalToNoisePeaks.OrderBy(peakSort => peakSort.X).ToList();
                                                    for (int i = 0; i < numIsotopologues; i++)
                                                    {
                                                        peaksFound[channelIndex + i, j] = secondTryOrderedByMzPeaks[i];
                                                    }
                                                }
                                                else
                                                {
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
                                        else
                                        {
                                            for (int i = 0; i < numIsotopologues; i++)
                                            {
                                                peaksFound[channelIndex + i, j] = null;
                                            }
                                        }
                                    }
                                    else
                                    {
                                        // Peaks found for an incomplete set of isotopologues by pattern detection
                                        if (topSignalToNoisePeaks.Count >= peaksNeeded)
                                        {
                                            ILabeledPeak[] mappedPeaks = mapPeaks(topSignalToNoisePeaks, channelIndex, channelIndex + (numIsotopologues - 1), j, charge);
                                            MassRange newSearchWindow = checkMappedPeaks(mappedPeaks, totalRange.Width, tolerance);

                                            if (newSearchWindow == null)
                                            {
                                                for (int m = channelIndex; m < channelIndex + numIsotopologues; m++)
                                                {
                                                    peaksFound[m, j] = mappedPeaks[m - channelIndex];
                                                }
                                            }
                                            else
                                            {
                                                List<ILabeledPeak> secondTryPeaks = findPeakPatterns(newSearchWindow, charge, current);
                                                List<ILabeledPeak> secondTryTopSignalToNoisePeaks = null;
                                                List<ILabeledPeak> secondTryOrderedByMzPeaks = null;
                                                if (secondTryPeaks != null && secondTryPeaks.Count > topSignalToNoisePeaks.Count)
                                                {
                                                    secondTryTopSignalToNoisePeaks = secondTryPeaks.OrderByDescending(peakSort => peakSort.GetSignalToNoise()).Take(numIsotopologues).ToList();
                                                    if (secondTryTopSignalToNoisePeaks.Count == numIsotopologues)
                                                    {
                                                        secondTryOrderedByMzPeaks = secondTryTopSignalToNoisePeaks.OrderBy(peakSort => peakSort.X).ToList();
                                                        for (int i = 0; i < numIsotopologues; i++)
                                                        {
                                                            peaksFound[channelIndex + i, j] = secondTryOrderedByMzPeaks[i];
                                                        }
                                                    }
                                                    else if (secondTryTopSignalToNoisePeaks.Count >= peaksNeeded)
                                                    {
                                                        double additionalPPMError = MassTolerance.GetTolerance(newSearchWindow.Minimum, totalRange.Minimum, MassToleranceType.PPM);
                                                        ILabeledPeak[] secondTryMappedPeaks = mapPeaks(secondTryTopSignalToNoisePeaks, channelIndex, channelIndex + (numIsotopologues - 1), j, charge, additionalPPMError);
                                                        for (int m = channelIndex; m < channelIndex + numIsotopologues; m++)
                                                        {
                                                            peaksFound[m, j] = secondTryMappedPeaks[m - channelIndex];
                                                        }
                                                    }
                                                }
                                                else
                                                {
                                                    for (int m = channelIndex; m < channelIndex + numIsotopologues; m++)
                                                    {
                                                        peaksFound[m, j] = mappedPeaks[m - channelIndex];
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }

                            // Count detected peaks
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
                                isotopeCount++;
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
                        if (isotopeCount > 0)
                        {
                            clusterCount++;
                        }
                    }
                } // End cluster loop

                if (clusterCount > 0)
                {
                    if (Form1.SEGMENTEDINJECTIONTIMES)
                    {
                        bool peakFound = false;
                        int channelCount = 0;
                        double mZ = 0;
                        double timeSum = 0;

                        while (channelCount < numChannels && !peakFound)
                        {
                            for (int j = 0; j < numIsotopes; j++)
                            {
                                if (peaksFound[channelCount, j] != null)
                                {
                                    mZ = peaksFound[channelCount, j].X;
                                    peakFound = true;
                                }
                            }                    
                            channelCount++;
                        }
                        
                        Dictionary<int, Dictionary<Range<double>, double>> scanNumbers = null;
                        if (Form1.INJECTIONTIMES.TryGetValue(rawFile.Name, out scanNumbers))
                        {
                            Dictionary<Range<double>, double> injectionTimes = null;

                            if (scanNumbers.TryGetValue(current.SpectrumNumber, out injectionTimes))
                            {
                                List<double> times = new List<double>();
                                foreach (KeyValuePair<Range<double>, double> segment in injectionTimes)
                                {
                                    if (segment.Key.Contains(mZ))
                                    {
                                        times.Add(segment.Value);
                                        timeSum += segment.Value;
                                    }
                                }

                                injectionTime = (timeSum / ((double)times.Count)) / 1000.0;
                            }
                        }
                    }
                    
                    pair = new Pair(this, rawFile, current.SpectrumNumber, injectionTime, current.RetentionTime);
                    pair.peaks = peaksFound;
                    if (pair.totalPeakCount > 0)
                    {
                        allPairs[rawFile].Add(pair);
                    }  
                }   
            }
            // Traditional SILAC (i.e., one channel per cluster)
            else
            {
                isotopeCount = 0;
                for (int j = 0; j < numIsotopes; j++)
                {
                    peaksCount = 0;
                    for (int i = 0; i < numChannels; i++)
                    {
                        peaksFound[i, j] = largestPeak(adjustedTheoMasses[i, j], current, TOLERANCE, rawFile);
                        if (peaksFound[i, j] != null)
                        {
                            peaksCount++;
                        }
                    }
                    if (peaksCount < peaksNeeded)
                    {
                        for (int k = 0; k < numChannels; k++)
                        {
                            peaksFound[k, j] = null;
                        }
                    }
                    else
                    {
                        isotopeCount++;
                    }
                }

                if (isotopeCount > 0)
                {
                    pair = new Pair(this, rawFile, current.SpectrumNumber, injectionTime, current.RetentionTime);
                    pair.peaks = peaksFound;
                    allPairs[rawFile].Add(pair);
                }
            }
        }

        public List<ILabeledPeak> cleanPeaks(List<ILabeledPeak> initialList, int charge)
        {
            if (initialList != null && initialList.Count > 1)
            {
                double theoTh = (spacingMassRange[1, 0].Mean * 1000.0) /  (double)charge;
                double theoThThreshold = theoTh / 1.5;
                List<ILabeledPeak> cleanedList = new List<ILabeledPeak>();
                List<ILabeledPeak> sortedIntensityList = initialList.OrderByDescending(dirtyPeak => dirtyPeak.GetSignalToNoise()).ToList();
                List<ILabeledPeak> potentialNubs = new List<ILabeledPeak>();

                for (int i = 0; i < sortedIntensityList.Count; i++)
                {
                    for (int j = i + 1; j < sortedIntensityList.Count; j++)
                    {
                        double mThDifference = Math.Abs(sortedIntensityList[j].X - sortedIntensityList[i].X) * 1000.0;
                        if (mThDifference < theoThThreshold)
                        {
                            if (!potentialNubs.Contains(sortedIntensityList[j])) potentialNubs.Add(sortedIntensityList[j]);
                        }
                    }
                }

                cleanedList = sortedIntensityList;
                foreach (ILabeledPeak peak in potentialNubs)
                {
                    cleanedList.Remove(peak);
                }

                return cleanedList;
            }
            else
            {
                return initialList;
            }
        }

        /* When given a list of peaks corresponding to an incomplete set of isotopologues, finds the combination with the smallest associated ppm error
         */
        public ILabeledPeak[] mapPeaks(List<ILabeledPeak> peaks, int channelStart, int channelEnd, int isotope, int charge, double additionalPPMError = 0.0)
        {
            // First perform any additional ppm adjustments
            double[,] additionallyAdjustedTheoMasses = adjustedTheoMasses;
            for (int i = channelStart; i <= channelEnd; i++)
            {
                additionallyAdjustedTheoMasses[i, isotope] = ((additionalPPMError * adjustedTheoMasses[i, isotope]) / 1000000.0) + adjustedTheoMasses[i, isotope];
            }
            
            int numChannelsToConsider = channelEnd - channelStart + 1;
            ILabeledPeak[] mappedPeaks = new ILabeledPeak[numChannelsToConsider];
            List<ILabeledPeak> orderedByMzPeaks = null;

            orderedByMzPeaks = peaks.OrderBy(peak => peak.X).ToList();

            int possibleCombinations = 0;
            if (orderedByMzPeaks.Count == 1)
            {
                if (numChannelsToConsider == 2) possibleCombinations = 2;
                else if (numChannelsToConsider == 3) possibleCombinations = 3;
            }
            else if (orderedByMzPeaks.Count == 2)
            {
                if (numChannelsToConsider == 3) possibleCombinations = 3;
                else if (numChannelsToConsider == 4) possibleCombinations = 6;
            }
            else if (orderedByMzPeaks.Count == 3)
            {
                if (numChannelsToConsider == 4) possibleCombinations = 4;
                else if (numChannelsToConsider == 6) possibleCombinations = 20;
            }
            else if (orderedByMzPeaks.Count == 4) possibleCombinations = 15;
            else if (orderedByMzPeaks.Count == 5) possibleCombinations = 6;
            KeyValuePair<double, ILabeledPeak>[,] errorCheck = new KeyValuePair<double, ILabeledPeak>[possibleCombinations, numChannelsToConsider];
            
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
            for (int i = 0; i < orderedByMzPeaks.Count; i++)
            {
                currentPeak = orderedByMzPeaks[i];

                if (numIsotopologues == 2)
                {
                    neutralMass = Mass.MassFromMz(currentPeak.X, charge);
                    error1 = MassTolerance.GetTolerance(neutralMass, additionallyAdjustedTheoMasses[channelStart, isotope], MassToleranceType.PPM);
                    error2 = MassTolerance.GetTolerance(neutralMass, additionallyAdjustedTheoMasses[channelEnd, isotope], MassToleranceType.PPM);
                    if ((error1) < (error2))
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
                    error1 = MassTolerance.GetTolerance(neutralMass, additionallyAdjustedTheoMasses[channelStart, isotope], MassToleranceType.PPM);
                    error2 = MassTolerance.GetTolerance(neutralMass, additionallyAdjustedTheoMasses[channelStart + 1, isotope], MassToleranceType.PPM);
                    error3 = MassTolerance.GetTolerance(neutralMass, additionallyAdjustedTheoMasses[channelStart + 2, isotope], MassToleranceType.PPM);

                    // 3 possibilities
                    if (orderedByMzPeaks.Count == 1)
                    {
                        if (error1 < error2 && error1 < error3) mappedPeaks[0] = currentPeak;
                        else if (error2 < error1 && error2 < error3) mappedPeaks[1] = currentPeak;
                        else if (error3 < error1 && error3 < error2) mappedPeaks[2] = currentPeak;
                    }
                    // 3 possibilities
                    else if (orderedByMzPeaks.Count == 2)
                    {
                        if (i == 0)
                        {
                            errorCheck[0, 0] = new KeyValuePair<double, ILabeledPeak>((error1), currentPeak);
                            errorCheck[1, 0] = new KeyValuePair<double, ILabeledPeak>((error1), currentPeak);
                            errorCheck[2, 1] = new KeyValuePair<double, ILabeledPeak>((error2), currentPeak);
                        }
                        else if (i == 1)
                        {
                            errorCheck[0, 1] = new KeyValuePair<double, ILabeledPeak>((error2), currentPeak);
                            errorCheck[1, 2] = new KeyValuePair<double, ILabeledPeak>((error3), currentPeak);
                            errorCheck[2, 2] = new KeyValuePair<double, ILabeledPeak>((error3), currentPeak);
                        }
                    }
                }
                else if (numIsotopologues == 4)
                {
                    neutralMass = Mass.MassFromMz(currentPeak.X, charge);
                    error1 = MassTolerance.GetTolerance(neutralMass, additionallyAdjustedTheoMasses[channelStart, isotope], MassToleranceType.PPM);
                    error2 = MassTolerance.GetTolerance(neutralMass, additionallyAdjustedTheoMasses[channelStart + 1, isotope], MassToleranceType.PPM);
                    error3 = MassTolerance.GetTolerance(neutralMass, additionallyAdjustedTheoMasses[channelStart + 2, isotope], MassToleranceType.PPM);
                    error4 = MassTolerance.GetTolerance(neutralMass, additionallyAdjustedTheoMasses[channelStart + 3, isotope], MassToleranceType.PPM);

                    // 6 possibilities
                    if (orderedByMzPeaks.Count == 2)
                    {
                        if (i == 0)
                        {
                            errorCheck[0, 0] = new KeyValuePair<double, ILabeledPeak>((error1), currentPeak);
                            errorCheck[1, 0] = new KeyValuePair<double, ILabeledPeak>((error1), currentPeak);
                            errorCheck[2, 0] = new KeyValuePair<double, ILabeledPeak>((error1), currentPeak);
                            errorCheck[3, 1] = new KeyValuePair<double, ILabeledPeak>((error2), currentPeak);
                            errorCheck[4, 1] = new KeyValuePair<double, ILabeledPeak>((error2), currentPeak);
                            errorCheck[5, 2] = new KeyValuePair<double, ILabeledPeak>((error3), currentPeak);
                        }
                        else
                        {
                            errorCheck[0, 1] = new KeyValuePair<double, ILabeledPeak>((error2), currentPeak);
                            errorCheck[1, 2] = new KeyValuePair<double, ILabeledPeak>((error3), currentPeak);
                            errorCheck[2, 3] = new KeyValuePair<double, ILabeledPeak>((error4), currentPeak);
                            errorCheck[3, 2] = new KeyValuePair<double, ILabeledPeak>((error3), currentPeak);
                            errorCheck[4, 3] = new KeyValuePair<double, ILabeledPeak>((error4), currentPeak);
                            errorCheck[5, 3] = new KeyValuePair<double, ILabeledPeak>((error4), currentPeak);
                        }
                    }

                    // 4 possibilities
                    if (orderedByMzPeaks.Count == 3)
                    {
                        if (i == 0)
                        {
                            errorCheck[0, 0] = new KeyValuePair<double, ILabeledPeak>((error1), currentPeak);
                            errorCheck[1, 0] = new KeyValuePair<double, ILabeledPeak>((error1), currentPeak);
                            errorCheck[2, 0] = new KeyValuePair<double, ILabeledPeak>((error1), currentPeak);
                            errorCheck[3, 1] = new KeyValuePair<double, ILabeledPeak>((error2), currentPeak);
                        }
                        else if (i == 1)
                        {
                            errorCheck[0, 1] = new KeyValuePair<double, ILabeledPeak>((error2), currentPeak);
                            errorCheck[1, 1] = new KeyValuePair<double, ILabeledPeak>((error2), currentPeak);
                            errorCheck[2, 2] = new KeyValuePair<double, ILabeledPeak>((error3), currentPeak);
                            errorCheck[3, 2] = new KeyValuePair<double, ILabeledPeak>((error3), currentPeak);
                        }
                        else
                        {
                            errorCheck[0, 2] = new KeyValuePair<double, ILabeledPeak>((error3), currentPeak);
                            errorCheck[1, 3] = new KeyValuePair<double, ILabeledPeak>((error4), currentPeak);
                            errorCheck[2, 3] = new KeyValuePair<double, ILabeledPeak>((error4), currentPeak);
                            errorCheck[3, 3] = new KeyValuePair<double, ILabeledPeak>((error4), currentPeak);
                        }
                    }
                }
                else if (numIsotopologues == 6)
                {
                    neutralMass = Mass.MassFromMz(currentPeak.X, charge);
                    error1 = MassTolerance.GetTolerance(neutralMass, additionallyAdjustedTheoMasses[channelStart, isotope], MassToleranceType.PPM);
                    error2 = MassTolerance.GetTolerance(neutralMass, additionallyAdjustedTheoMasses[channelStart + 1, isotope], MassToleranceType.PPM);
                    error3 = MassTolerance.GetTolerance(neutralMass, additionallyAdjustedTheoMasses[channelStart + 2, isotope], MassToleranceType.PPM);
                    error4 = MassTolerance.GetTolerance(neutralMass, additionallyAdjustedTheoMasses[channelStart + 3, isotope], MassToleranceType.PPM);
                    error5 = MassTolerance.GetTolerance(neutralMass, additionallyAdjustedTheoMasses[channelStart + 4, isotope], MassToleranceType.PPM);
                    error6 = MassTolerance.GetTolerance(neutralMass, additionallyAdjustedTheoMasses[channelStart + 5, isotope], MassToleranceType.PPM);

                    // 20 possibilities
                    if (orderedByMzPeaks.Count == 3)
                    {
                        if (i == 0)
                        {
                            errorCheck[0, 0] = new KeyValuePair<double, ILabeledPeak>((error1), currentPeak);
                            errorCheck[1, 0] = new KeyValuePair<double, ILabeledPeak>((error1), currentPeak);
                            errorCheck[2, 0] = new KeyValuePair<double, ILabeledPeak>((error1), currentPeak);
                            errorCheck[3, 0] = new KeyValuePair<double, ILabeledPeak>((error1), currentPeak);
                            errorCheck[4, 0] = new KeyValuePair<double, ILabeledPeak>((error1), currentPeak);
                            errorCheck[5, 0] = new KeyValuePair<double, ILabeledPeak>((error1), currentPeak);
                            errorCheck[6, 0] = new KeyValuePair<double, ILabeledPeak>((error1), currentPeak);
                            errorCheck[7, 0] = new KeyValuePair<double, ILabeledPeak>((error1), currentPeak);
                            errorCheck[8, 0] = new KeyValuePair<double, ILabeledPeak>((error1), currentPeak);
                            errorCheck[9, 0] = new KeyValuePair<double, ILabeledPeak>((error1), currentPeak);
                            errorCheck[10, 1] = new KeyValuePair<double, ILabeledPeak>((error2), currentPeak);
                            errorCheck[11, 1] = new KeyValuePair<double, ILabeledPeak>((error2), currentPeak);
                            errorCheck[12, 1] = new KeyValuePair<double, ILabeledPeak>((error2), currentPeak);
                            errorCheck[13, 1] = new KeyValuePair<double, ILabeledPeak>((error2), currentPeak);
                            errorCheck[14, 1] = new KeyValuePair<double, ILabeledPeak>((error2), currentPeak);
                            errorCheck[15, 1] = new KeyValuePair<double, ILabeledPeak>((error2), currentPeak);
                            errorCheck[16, 2] = new KeyValuePair<double, ILabeledPeak>((error3), currentPeak);
                            errorCheck[17, 2] = new KeyValuePair<double, ILabeledPeak>((error3), currentPeak);
                            errorCheck[18, 2] = new KeyValuePair<double, ILabeledPeak>((error3), currentPeak);
                            errorCheck[19, 3] = new KeyValuePair<double, ILabeledPeak>((error4), currentPeak);
                            
                        }
                        else if (i == 1)
                        {
                            errorCheck[0, 1] = new KeyValuePair<double, ILabeledPeak>((error2), currentPeak);
                            errorCheck[1, 1] = new KeyValuePair<double, ILabeledPeak>((error2), currentPeak);
                            errorCheck[2, 1] = new KeyValuePair<double, ILabeledPeak>((error2), currentPeak);
                            errorCheck[3, 1] = new KeyValuePair<double, ILabeledPeak>((error2), currentPeak);
                            errorCheck[4, 2] = new KeyValuePair<double, ILabeledPeak>((error3), currentPeak);
                            errorCheck[5, 2] = new KeyValuePair<double, ILabeledPeak>((error3), currentPeak);
                            errorCheck[6, 2] = new KeyValuePair<double, ILabeledPeak>((error3), currentPeak);
                            errorCheck[7, 3] = new KeyValuePair<double, ILabeledPeak>((error4), currentPeak);
                            errorCheck[8, 3] = new KeyValuePair<double, ILabeledPeak>((error4), currentPeak);
                            errorCheck[9, 4] = new KeyValuePair<double, ILabeledPeak>((error5), currentPeak);
                            errorCheck[10, 2] = new KeyValuePair<double, ILabeledPeak>((error3), currentPeak);
                            errorCheck[11, 2] = new KeyValuePair<double, ILabeledPeak>((error3), currentPeak);
                            errorCheck[12, 2] = new KeyValuePair<double, ILabeledPeak>((error3), currentPeak);
                            errorCheck[13, 3] = new KeyValuePair<double, ILabeledPeak>((error4), currentPeak);
                            errorCheck[14, 3] = new KeyValuePair<double, ILabeledPeak>((error4), currentPeak);
                            errorCheck[15, 4] = new KeyValuePair<double, ILabeledPeak>((error5), currentPeak);
                            errorCheck[16, 3] = new KeyValuePair<double, ILabeledPeak>((error4), currentPeak);
                            errorCheck[17, 3] = new KeyValuePair<double, ILabeledPeak>((error4), currentPeak);
                            errorCheck[18, 4] = new KeyValuePair<double, ILabeledPeak>((error5), currentPeak);
                            errorCheck[19, 4] = new KeyValuePair<double, ILabeledPeak>((error5), currentPeak);
                        }
                        else if (i == 2)
                        {
                            errorCheck[0, 2] = new KeyValuePair<double, ILabeledPeak>((error3), currentPeak);
                            errorCheck[1, 3] = new KeyValuePair<double, ILabeledPeak>((error4), currentPeak);
                            errorCheck[2, 4] = new KeyValuePair<double, ILabeledPeak>((error5), currentPeak);
                            errorCheck[3, 5] = new KeyValuePair<double, ILabeledPeak>((error6), currentPeak);
                            errorCheck[4, 3] = new KeyValuePair<double, ILabeledPeak>((error4), currentPeak);
                            errorCheck[5, 4] = new KeyValuePair<double, ILabeledPeak>((error5), currentPeak);
                            errorCheck[6, 5] = new KeyValuePair<double, ILabeledPeak>((error6), currentPeak);
                            errorCheck[7, 4] = new KeyValuePair<double, ILabeledPeak>((error5), currentPeak);
                            errorCheck[8, 5] = new KeyValuePair<double, ILabeledPeak>((error6), currentPeak);
                            errorCheck[9, 5] = new KeyValuePair<double, ILabeledPeak>((error6), currentPeak);
                            errorCheck[10, 3] = new KeyValuePair<double, ILabeledPeak>((error4), currentPeak);
                            errorCheck[11, 4] = new KeyValuePair<double, ILabeledPeak>((error5), currentPeak);
                            errorCheck[12, 5] = new KeyValuePair<double, ILabeledPeak>((error6), currentPeak);
                            errorCheck[13, 4] = new KeyValuePair<double, ILabeledPeak>((error5), currentPeak);
                            errorCheck[14, 5] = new KeyValuePair<double, ILabeledPeak>((error6), currentPeak);
                            errorCheck[15, 5] = new KeyValuePair<double, ILabeledPeak>((error6), currentPeak);
                            errorCheck[16, 4] = new KeyValuePair<double, ILabeledPeak>((error5), currentPeak);
                            errorCheck[17, 5] = new KeyValuePair<double, ILabeledPeak>((error6), currentPeak);
                            errorCheck[18, 5] = new KeyValuePair<double, ILabeledPeak>((error6), currentPeak);
                            errorCheck[19, 5] = new KeyValuePair<double, ILabeledPeak>((error6), currentPeak);
                        }
                    }
                    // 15 possibilities
                    else if (orderedByMzPeaks.Count == 4)
                    {
                        if (i == 0)
                        {
                            errorCheck[0, 0] = new KeyValuePair<double, ILabeledPeak>((error1), currentPeak);
                            errorCheck[1, 0] = new KeyValuePair<double, ILabeledPeak>((error1), currentPeak);
                            errorCheck[2, 0] = new KeyValuePair<double, ILabeledPeak>((error1), currentPeak);
                            errorCheck[3, 0] = new KeyValuePair<double, ILabeledPeak>((error1), currentPeak);
                            errorCheck[4, 0] = new KeyValuePair<double, ILabeledPeak>((error1), currentPeak);
                            errorCheck[5, 0] = new KeyValuePair<double, ILabeledPeak>((error2), currentPeak);
                            errorCheck[6, 0] = new KeyValuePair<double, ILabeledPeak>((error1), currentPeak);
                            errorCheck[7, 0] = new KeyValuePair<double, ILabeledPeak>((error1), currentPeak);
                            errorCheck[8, 0] = new KeyValuePair<double, ILabeledPeak>((error1), currentPeak);
                            errorCheck[9, 0] = new KeyValuePair<double, ILabeledPeak>((error1), currentPeak);
                            errorCheck[10, 1] = new KeyValuePair<double, ILabeledPeak>((error2), currentPeak);
                            errorCheck[11, 1] = new KeyValuePair<double, ILabeledPeak>((error2), currentPeak);
                            errorCheck[12, 1] = new KeyValuePair<double, ILabeledPeak>((error2), currentPeak);
                            errorCheck[13, 1] = new KeyValuePair<double, ILabeledPeak>((error2), currentPeak);
                            errorCheck[14, 2] = new KeyValuePair<double, ILabeledPeak>((error3), currentPeak);                            
                        }
                        else if (i == 1)
                        {
                            errorCheck[0, 1] = new KeyValuePair<double, ILabeledPeak>((error2), currentPeak);
                            errorCheck[1, 1] = new KeyValuePair<double, ILabeledPeak>((error2), currentPeak);
                            errorCheck[2, 1] = new KeyValuePair<double, ILabeledPeak>((error2), currentPeak);
                            errorCheck[3, 1] = new KeyValuePair<double, ILabeledPeak>((error2), currentPeak);
                            errorCheck[4, 1] = new KeyValuePair<double, ILabeledPeak>((error2), currentPeak);
                            errorCheck[5, 1] = new KeyValuePair<double, ILabeledPeak>((error2), currentPeak);
                            errorCheck[6, 2] = new KeyValuePair<double, ILabeledPeak>((error3), currentPeak);
                            errorCheck[7, 2] = new KeyValuePair<double, ILabeledPeak>((error3), currentPeak);
                            errorCheck[8, 2] = new KeyValuePair<double, ILabeledPeak>((error3), currentPeak);
                            errorCheck[9, 3] = new KeyValuePair<double, ILabeledPeak>((error4), currentPeak);
                            errorCheck[10, 2] = new KeyValuePair<double, ILabeledPeak>((error3), currentPeak);
                            errorCheck[11, 2] = new KeyValuePair<double, ILabeledPeak>((error3), currentPeak);
                            errorCheck[12, 2] = new KeyValuePair<double, ILabeledPeak>((error3), currentPeak);
                            errorCheck[13, 3] = new KeyValuePair<double, ILabeledPeak>((error4), currentPeak);
                            errorCheck[14, 3] = new KeyValuePair<double, ILabeledPeak>((error4), currentPeak);
                        }
                        else if (i == 2)
                        {
                            errorCheck[0, 2] = new KeyValuePair<double, ILabeledPeak>((error3), currentPeak);
                            errorCheck[1, 2] = new KeyValuePair<double, ILabeledPeak>((error3), currentPeak);
                            errorCheck[2, 2] = new KeyValuePair<double, ILabeledPeak>((error3), currentPeak);
                            errorCheck[3, 3] = new KeyValuePair<double, ILabeledPeak>((error4), currentPeak);
                            errorCheck[4, 3] = new KeyValuePair<double, ILabeledPeak>((error4), currentPeak);
                            errorCheck[5, 4] = new KeyValuePair<double, ILabeledPeak>((error5), currentPeak);
                            errorCheck[6, 3] = new KeyValuePair<double, ILabeledPeak>((error4), currentPeak);
                            errorCheck[7, 3] = new KeyValuePair<double, ILabeledPeak>((error4), currentPeak);
                            errorCheck[8, 4] = new KeyValuePair<double, ILabeledPeak>((error5), currentPeak);
                            errorCheck[9, 4] = new KeyValuePair<double, ILabeledPeak>((error5), currentPeak);
                            errorCheck[10, 3] = new KeyValuePair<double, ILabeledPeak>((error4), currentPeak);
                            errorCheck[11, 3] = new KeyValuePair<double, ILabeledPeak>((error4), currentPeak);
                            errorCheck[12, 4] = new KeyValuePair<double, ILabeledPeak>((error5), currentPeak);
                            errorCheck[13, 4] = new KeyValuePair<double, ILabeledPeak>((error5), currentPeak);
                            errorCheck[14, 4] = new KeyValuePair<double, ILabeledPeak>((error5), currentPeak);
                        }
                        else if (i == 3)
                        {
                            errorCheck[0, 3] = new KeyValuePair<double, ILabeledPeak>((error4), currentPeak);
                            errorCheck[1, 4] = new KeyValuePair<double, ILabeledPeak>((error5), currentPeak);
                            errorCheck[2, 5] = new KeyValuePair<double, ILabeledPeak>((error6), currentPeak);
                            errorCheck[3, 4] = new KeyValuePair<double, ILabeledPeak>((error5), currentPeak);
                            errorCheck[4, 5] = new KeyValuePair<double, ILabeledPeak>((error6), currentPeak);
                            errorCheck[5, 5] = new KeyValuePair<double, ILabeledPeak>((error6), currentPeak);
                            errorCheck[6, 4] = new KeyValuePair<double, ILabeledPeak>((error5), currentPeak);
                            errorCheck[7, 5] = new KeyValuePair<double, ILabeledPeak>((error6), currentPeak);
                            errorCheck[8, 5] = new KeyValuePair<double, ILabeledPeak>((error6), currentPeak);
                            errorCheck[9, 5] = new KeyValuePair<double, ILabeledPeak>((error6), currentPeak);
                            errorCheck[10, 4] = new KeyValuePair<double, ILabeledPeak>((error5), currentPeak);
                            errorCheck[11, 5] = new KeyValuePair<double, ILabeledPeak>((error6), currentPeak);
                            errorCheck[12, 5] = new KeyValuePair<double, ILabeledPeak>((error6), currentPeak);
                            errorCheck[13, 5] = new KeyValuePair<double, ILabeledPeak>((error6), currentPeak);
                            errorCheck[14, 5] = new KeyValuePair<double, ILabeledPeak>((error6), currentPeak);
                        }
                    }
                    // 6 possibilities
                    else if (orderedByMzPeaks.Count == 5)
                    {
                        if (i == 0)
                        {
                            errorCheck[0, 0] = new KeyValuePair<double, ILabeledPeak>((error1), currentPeak);
                            errorCheck[1, 0] = new KeyValuePair<double, ILabeledPeak>((error1), currentPeak);
                            errorCheck[2, 0] = new KeyValuePair<double, ILabeledPeak>((error1), currentPeak);
                            errorCheck[3, 0] = new KeyValuePair<double, ILabeledPeak>((error1), currentPeak);
                            errorCheck[4, 0] = new KeyValuePair<double, ILabeledPeak>((error1), currentPeak);
                            errorCheck[5, 1] = new KeyValuePair<double, ILabeledPeak>((error2), currentPeak);
                        }
                        else if (i == 1)
                        {
                            errorCheck[0, 1] = new KeyValuePair<double, ILabeledPeak>((error2), currentPeak);
                            errorCheck[1, 1] = new KeyValuePair<double, ILabeledPeak>((error2), currentPeak);
                            errorCheck[2, 1] = new KeyValuePair<double, ILabeledPeak>((error2), currentPeak);
                            errorCheck[3, 1] = new KeyValuePair<double, ILabeledPeak>((error2), currentPeak);
                            errorCheck[4, 2] = new KeyValuePair<double, ILabeledPeak>((error3), currentPeak);
                            errorCheck[5, 2] = new KeyValuePair<double, ILabeledPeak>((error3), currentPeak);
                        }
                        else if (i == 2)
                        {
                            errorCheck[0, 2] = new KeyValuePair<double, ILabeledPeak>((error3), currentPeak);
                            errorCheck[1, 2] = new KeyValuePair<double, ILabeledPeak>((error3), currentPeak);
                            errorCheck[2, 2] = new KeyValuePair<double, ILabeledPeak>((error3), currentPeak);
                            errorCheck[3, 3] = new KeyValuePair<double, ILabeledPeak>((error4), currentPeak);
                            errorCheck[4, 3] = new KeyValuePair<double, ILabeledPeak>((error4), currentPeak);
                            errorCheck[5, 3] = new KeyValuePair<double, ILabeledPeak>((error4), currentPeak);
                        }
                        else if (i == 3)
                        {
                            errorCheck[0, 3] = new KeyValuePair<double, ILabeledPeak>((error4), currentPeak);
                            errorCheck[1, 3] = new KeyValuePair<double, ILabeledPeak>((error4), currentPeak);
                            errorCheck[2, 4] = new KeyValuePair<double, ILabeledPeak>((error5), currentPeak);
                            errorCheck[3, 4] = new KeyValuePair<double, ILabeledPeak>((error5), currentPeak);
                            errorCheck[4, 4] = new KeyValuePair<double, ILabeledPeak>((error5), currentPeak);
                            errorCheck[5, 4] = new KeyValuePair<double, ILabeledPeak>((error5), currentPeak);
                        }
                        else if (i == 4)
                        {
                            errorCheck[0, 4] = new KeyValuePair<double, ILabeledPeak>((error5), currentPeak);
                            errorCheck[1, 5] = new KeyValuePair<double, ILabeledPeak>((error6), currentPeak);
                            errorCheck[2, 5] = new KeyValuePair<double, ILabeledPeak>((error6), currentPeak);
                            errorCheck[3, 5] = new KeyValuePair<double, ILabeledPeak>((error6), currentPeak);
                            errorCheck[4, 5] = new KeyValuePair<double, ILabeledPeak>((error6), currentPeak);
                            errorCheck[5, 5] = new KeyValuePair<double, ILabeledPeak>((error6), currentPeak);
                        }
                    }
                }
            }

            // Find combination that leads to lowest overall error
            List<double> individualPPMErrors = new List<double>();
            double errorRange;
            minError = 1000;
            minErrorIndex = -1;
            for (int j = 0; j < possibleCombinations; j++)
            {
                overallError = 0;
                for (int k = 0; k < numChannelsToConsider; k++)
                {
                    overallError += Math.Abs(errorCheck[j, k].Key);
                    individualPPMErrors.Add(errorCheck[j, k].Key);
                }

                individualPPMErrors.Sort();
                errorRange = Math.Abs(individualPPMErrors[individualPPMErrors.Count - 1] - individualPPMErrors[0]);
                overallError += errorRange;

                if (overallError < minError)
                {
                    minError = overallError;
                    minErrorIndex = j;
                }
            }

            // Set mappedPeaks to peak combination with lowest overall error
            for (int n = 0; n < numChannelsToConsider; n++)
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
                            //else if (log10MaxIntensity <= System.Math.Log10(Form1.MAXIMUMDNL))
                            //{
                                
                            //}
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
                            //else if (log10MaxIntensity <= System.Math.Log10(Form1.MAXIMUMDNL))
                            //{
                                
                            //}
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
        public void checkPairSpacing(MSDataFile rawFile, List<Spacing> spacings = null)
        {            
            int charge = bestPSMs[rawFile].Charge;
            List<Pair> pairs;
            allPairs.TryGetValue(rawFile, out pairs);
            Pair pair;
            int peakCount;
            bool[,] spacingChecks;
            bool spacingChecked;

            if (numIsotopologues < 2 && numClusterLabels > 0)
            {
                for (int i = 0; i < pairs.Count; i++)
                {
                    pair = pairs[i];

                    for (int k = 0; k < numIsotopes; k++)
                    {
                        peakCount = pair.peakCount[0, k]; ;
                        spacingChecks = pair.GetSpacingCheckArray(0, k);

                        for (int j = 0; j < numChannels; j++)
                        {
                            spacingChecked = false;
                            for (int c = 0; c < numChannels; c++)
                            {
                                if (spacingChecks[j, c]) spacingChecked = true;
                            }
                            if (pair.peaks[j, k] != null && !spacingChecked)
                            {
                                peakCount--;
                                pair.peaks[j, k] = null;
                            }
                        }

                        if (peakCount < peaksNeeded)
                        {
                            for (int j = 0; j < numChannels; j++)
                            {
                                pair.peaks[j, k] = null;
                            }
                        }
                    }
                }
            }
            else if (numIsotopologues > 1 && numIsotopologueLabels > 0)
            {
                for (int i = 0; i < pairs.Count; i++)
                {
                    pair = pairs[i];

                    for (int c = 0; c < numClusters; c++)
                    {
                        for (int k = 0; k < numIsotopes; k++)
                        {
                            peakCount = pair.peakCount[c, k];
                            spacingChecks = pair.GetSpacingCheckArray(c, k);
                            int channelIndex = c * numIsotopologues;

                            for (int j = 0; j < numIsotopologues; j++)
                            {
                                spacingChecked = false;
                                for (int m = 0; m < numIsotopologues; m++)
                                {
                                    if (spacingChecks[j, m]) spacingChecked = true;
                                }
                                if (pair.peaks[j + channelIndex, k] != null && !spacingChecked)
                                {
                                    peakCount--;
                                    pair.peaks[j + channelIndex, k] = null;
                                }
                            }

                            if (peakCount < peaksNeeded)
                            {
                                for (int j = channelIndex; j < channelIndex + numIsotopologues; j++)
                                {
                                    pair.peaks[j, k] = null;
                                }
                            }
                        }
                    }
                }
            }

                    //if (numChannels == 2)
            //{
            //    for (int k = 0; k < numIsotopes; k++)
            //    {
            //        peakCount = 0;
            //        peak1 = pair.peaks[0, k];
            //        peak2 = pair.peaks[1, k];

                    //        if (peak1 != null) peakCount++;
            //        if (peak2 != null) peakCount++;

                    //        if (peakCount == 2)
            //        {
            //            peak1NeutralMass = Mass.MassFromMz(peak1.X, charge);
            //            peak2NeutralMass = Mass.MassFromMz(peak2.X, charge);

                    //            spacing1 = peak2NeutralMass - peak1NeutralMass;

                    //            //Remove pairs whose spacing is more than 10 mDa away from the calculated spacing
            //            if (spacing1 < spacingMassRange[0].Minimum || spacing1 > spacingMassRange[0].Maximum)
            //            {
            //                //Console.WriteLine("bad spacing");
            //                pair.peaks[0, k] = null;
            //                peakCount--;
            //                pair.peaks[1, k] = null;
            //                peakCount--;
            //            }
            //        }
            //    }
            //}
            //else if (numChannels == 3)
            //{
            //    for (int k = 0; k < numIsotopes; k++)
            //    {
            //        peakCount = 0;
            //        peak1 = pair.peaks[0, k];
            //        peak2 = pair.peaks[1, k];
            //        peak3 = pair.peaks[2, k];

                    //        if (peak1 != null) peakCount++;
            //        if (peak2 != null) peakCount++;
            //        if (peak3 != null) peakCount++;

                    //        peakCount = 0;
            //        peak1 = pair.peaks[0, k];
            //        peak2 = pair.peaks[1, k];
            //        peak3 = pair.peaks[2, k];

                    //        if (peak1 != null) peakCount++;
            //        if (peak2 != null) peakCount++;
            //        if (peak3 != null) peakCount++;

                    //        if (peakCount == 3)
            //        {
            //            peak1NeutralMass = Mass.MassFromMz(peak1.X, charge);
            //            peak2NeutralMass = Mass.MassFromMz(peak2.X, charge);
            //            peak3NeutralMass = Mass.MassFromMz(peak3.X, charge);

                    //            spacing1 = peak2NeutralMass - peak1NeutralMass;
            //            spacing2 = peak3NeutralMass - peak2NeutralMass;

                    //            if (Form1.NOISEBANDCAP)
            //            {
            //                if (!spacingMassRange[0].Contains(spacing1))
            //                {
            //                    pair.peaks[0, k] = null;
            //                    peakCount--;

                    //                    if (!spacingMassRange[1].Contains(spacing2))
            //                    {
            //                        pair.peaks[1, k] = null;
            //                        peakCount--;
            //                        pair.peaks[2, k] = null;
            //                        peakCount--;
            //                    }
            //                }
            //                else if (!spacingMassRange[1].Contains(spacing2))
            //                {
            //                    pair.peaks[2, k] = null;
            //                    peakCount--;
            //                }

                    //                if (peakCount < peaksNeeded)
            //                {
            //                    pair.peaks[0, k] = null;
            //                    pair.peaks[1, k] = null;
            //                    pair.peaks[2, k] = null;
            //                }
            //            }
            //            else
            //            {
            //                if (!spacingMassRange[0].Contains(spacing1) || !spacingMassRange[1].Contains(spacing2))
            //                {
            //                    pair.peaks[0, k] = null;
            //                    pair.peaks[1, k] = null;
            //                    pair.peaks[2, k] = null;
            //                }
            //            }
            //        }
            //        else if (Form1.NOISEBANDCAP && peakCount == 2)
            //        {
            //            if (peak1 != null)
            //            {
            //                peak1NeutralMass = Mass.MassFromMz(peak1.X, charge);

                    //                // Peaks 1 & 2
            //                if (peak2 != null)
            //                {
            //                    peak2NeutralMass = Mass.MassFromMz(peak2.X, charge);
            //                    spacing1 = peak2NeutralMass - peak1NeutralMass;
            //                    if (!spacingMassRange[0].Contains(spacing1))
            //                    {
            //                        pair.peaks[0, k] = null;
            //                        pair.peaks[1, k] = null;
            //                    }
            //                }
            //                // Peaks 1 & 3
            //                else if (peak3 != null)
            //                {
            //                    peak2NeutralMass = Mass.MassFromMz(peak3.X, charge);
            //                    spacing1 = peak2NeutralMass - peak1NeutralMass;
            //                    double extendedTheoSpacing = spacingMassRange[0].Mean + spacingMassRange[1].Mean;
            //                    MassRange doubleSpacing = new MassRange(extendedTheoSpacing, new MassTolerance(MassToleranceType.DA, SPACINGPERCENTAGEERROR * extendedTheoSpacing));
            //                    if (!doubleSpacing.Contains(spacing1))
            //                    {
            //                        pair.peaks[0, k] = null;
            //                        pair.peaks[2, k] = null;
            //                    }
            //                }
            //            }
            //            // Peaks 2 & 3
            //            else if (peak2 != null)
            //            {
            //                peak1NeutralMass = Mass.MassFromMz(peak2.X, charge);
            //                peak2NeutralMass = Mass.MassFromMz(peak3.X, charge);
            //                spacing1 = peak2NeutralMass - peak1NeutralMass;
            //                if (!spacingMassRange[1].Contains(spacing1))
            //                {
            //                    pair.peaks[1, k] = null;
            //                    pair.peaks[2, k] = null;
            //                }
            //            }
            //        }
            //        else if (Form1.NOISEBANDCAP && peakCount < peaksNeeded)
            //        {
            //            pair.peaks[0, k] = null;
            //            pair.peaks[1, k] = null;
            //            pair.peaks[2, k] = null;
            //        }

                    //        if (peakCount < peaksNeeded)
            //        {
            //            pair.peaks[0, k] = null;
            //            pair.peaks[1, k] = null;
            //            pair.peaks[2, k] = null;
            //        }
            //    }
            //}
            //else
            //{
            //    for (int i = 0; i < pairs.Count; i++)
            //    {
            //        pair = pairs[i];

            //        for (int k = 0; k < numIsotopes; k++)
            //        {
            //            peakCount = 0;
            //            peakSpacings = new double[numIsotopologues - 1];

            //            for (int c = 0; c < numClusters; c++)
            //            {
            //                int channelIndex = c * numIsotopologues;

            //                // Count # of peaks
            //                for (int j = channelIndex; j < channelIndex + numIsotopologues; j++)
            //                {
            //                    if (pair.peaks[j, k] != null) peakCount++;
            //                }

            //                // Complete peak sets
            //                if (peakCount == numChannels)
            //                {
            //                    for (int j = channelIndex; j < channelIndex + numIsotopologues - 1; j++)
            //                    {
            //                        peakSpacings[j] = Mass.MassFromMz(pair.peaks[j + 1, k].X - pair.peaks[j, k].X, charge);
            //                    }

            //                    if (peakSpacings.Count() > 1)
            //                    {
            //                        outerPeakSpacing = peakSpacings[peakSpacings.Count() - 1] + peakSpacings[0];
            //                        outerPeakSpacingRange = new MassRange(outerPeakSpacing, new MassTolerance(MassToleranceType.DA, SPACINGPERCENTAGEERROR * outerPeakSpacing));
            //                    }

            //                    if (numIsotopologues == 2)
            //                    {
            //                        if (!spacingMassRange[0].Contains(peakSpacings[0]))
            //                        {
            //                            pair.peaks[channelIndex, k] = null;
            //                            pair.peaks[channelIndex + 1, k] = null;
            //                            peakCount -= 2;
            //                        }
            //                    }
            //                    else if (numIsotopologues == 3)
            //                    {
            //                        // Bad outer spacing
            //                        if (!outerPeakSpacingRange.Contains(outerPeakSpacing))
            //                        {
            //                            // First inner spacing bad
            //                            if (!spacingMassRange[0].Contains(peakSpacings[0]))
            //                            {
            //                                pair.peaks[channelIndex, k] = null;
            //                                peakCount--;

            //                                // Both inner spacings bad
            //                                if (!spacingMassRange[1].Contains(peakSpacings[1]))
            //                                {
            //                                    pair.peaks[channelIndex + 1, k] = null;
            //                                    pair.peaks[channelIndex + 2, k] = null;
            //                                    peakCount -= 2;
            //                                }
            //                            }
            //                            // Second inner spacing bad
            //                            else if (!spacingMassRange[1].Contains(peakSpacings[1]))
            //                            {
            //                                pair.peaks[channelIndex + 2, k] = null;
            //                                peakCount--;
            //                            }
            //                        }
            //                        // Good outer spacing
            //                        else
            //                        {
            //                            // Both inner spacings bad
            //                            if (!spacingMassRange[0].Contains(peakSpacings[0]) && !spacingMassRange[1].Contains(peakSpacings[1]))
            //                            {
            //                                pair.peaks[channelIndex + 1, k] = null;
            //                                peakCount--;
            //                            }
            //                        }
            //                    }
            //                    else if (numIsotopologues == 4)
            //                    {
            //                        // Bad outer spacing
            //                        if (!outerPeakSpacingRange.Contains(outerPeakSpacing))
            //                        {
            //                            // First inner spacing bad
            //                            if (!spacingMassRange[0].Contains(peakSpacings[0]))
            //                            {
            //                                pair.peaks[channelIndex, k] = null;
            //                                peakCount--;

            //                                // First & second inner spacings bad
            //                                if (!spacingMassRange[1].Contains(peakSpacings[1]))
            //                                {
            //                                    pair.peaks[channelIndex + 1, k] = null;
            //                                    peakCount--;

            //                                    // All inner spacings bad
            //                                    if (!spacingMassRange[2].Contains(peakSpacings[2]))
            //                                    {
            //                                        pair.peaks[channelIndex + 2, k] = null;
            //                                        pair.peaks[channelIndex + 3, k] = null;
            //                                        peakCount -= 2;
            //                                    }
            //                                }

            //                                // First & third inner spacings bad
            //                                if (!spacingMassRange[2].Contains(peakSpacings[2]))
            //                                {
            //                                    pair.peaks[channelIndex + 3, k] = null;
            //                                    peakCount--;
            //                                }
            //                            }
            //                            // Second inner spacing bad
            //                            else if (!spacingMassRange[1].Contains(peakSpacings[1]))
            //                            {
            //                                // Second & third inner spacings bad
            //                                if (!spacingMassRange[2].Contains(peakSpacings[2]))
            //                                {
            //                                    pair.peaks[channelIndex + 2, k] = null;
            //                                    pair.peaks[channelIndex + 3, k] = null;
            //                                    peakCount -= 2;
            //                                }
            //                            }
            //                            // Third inner spacing bad
            //                            else if (!spacingMassRange[2].Contains(peakSpacings[2]))
            //                            {
            //                                pair.peaks[channelIndex + 3, k] = null;
            //                                peakCount--;
            //                            }
            //                        }
            //                        // Good outer spacing
            //                        else
            //                        {
            //                            // First inner spacing bad
            //                            if (!spacingMassRange[0].Contains(peakSpacings[0]))
            //                            {
            //                                // First & second inner spacings bad
            //                                if (!spacingMassRange[1].Contains(peakSpacings[1]))
            //                                {
            //                                    pair.peaks[channelIndex + 1, k] = null;
            //                                    peakCount--;

            //                                    // All inner spacings bad
            //                                    if (!spacingMassRange[2].Contains(peakSpacings[2]))
            //                                    {
            //                                        pair.peaks[channelIndex + 2, k] = null;
            //                                        peakCount--;
            //                                    }
            //                                }
            //                                // First & third inner spacings bad
            //                                if (!spacingMassRange[2].Contains(peakSpacings[2])) { }
            //                            }
            //                            // Second inner spacing bad
            //                            else if (!spacingMassRange[1].Contains(peakSpacings[1]))
            //                            {
            //                                // Second & third inner spacings bad
            //                                if (!spacingMassRange[2].Contains(peakSpacings[2]))
            //                                {
            //                                    pair.peaks[channelIndex + 2, k] = null;
            //                                    peakCount--;
            //                                }
            //                            }
            //                            // Third inner spacing bad
            //                            else if (!spacingMassRange[2].Contains(peakSpacings[2])) { }
            //                        }
            //                    }
            //                    else if (numIsotopologues == 6)
            //                    {
            //                        // Bad outer spacing
            //                        if (!outerPeakSpacingRange.Contains(outerPeakSpacing))
            //                        {
            //                            // First inner spacing bad
            //                            if (!spacingMassRange[0].Contains(peakSpacings[0]))
            //                            {
            //                                pair.peaks[channelIndex, k] = null;
            //                                peakCount--;

            //                                // First & second inner spacings bad
            //                                if (!spacingMassRange[1].Contains(peakSpacings[1]))
            //                                {
            //                                    pair.peaks[channelIndex + 1, k] = null;
            //                                    peakCount--;

            //                                    // All inner spacings bad
            //                                    if (!spacingMassRange[2].Contains(peakSpacings[2]))
            //                                    {
            //                                        pair.peaks[channelIndex + 2, k] = null;
            //                                        pair.peaks[channelIndex + 3, k] = null;
            //                                        peakCount -= 2;
            //                                    }
            //                                }

            //                                // First & third inner spacings bad
            //                                if (!spacingMassRange[2].Contains(peakSpacings[2]))
            //                                {
            //                                    pair.peaks[channelIndex + 3, k] = null;
            //                                    peakCount--;
            //                                }
            //                            }
            //                            // Second inner spacing bad
            //                            else if (!spacingMassRange[1].Contains(peakSpacings[1]))
            //                            {
            //                                // Second & third inner spacings bad
            //                                if (!spacingMassRange[2].Contains(peakSpacings[2]))
            //                                {
            //                                    pair.peaks[channelIndex + 2, k] = null;
            //                                    pair.peaks[channelIndex + 3, k] = null;
            //                                    peakCount -= 2;
            //                                }
            //                            }
            //                            // Third inner spacing bad
            //                            else if (!spacingMassRange[2].Contains(peakSpacings[2]))
            //                            {
            //                                pair.peaks[channelIndex + 3, k] = null;
            //                                peakCount--;
            //                            }
            //                        }
            //                        // Good outer spacing
            //                        else
            //                        {
            //                            // First inner spacing bad
            //                            if (!spacingMassRange[0].Contains(peakSpacings[0]))
            //                            {
            //                                // First & second inner spacings bad
            //                                if (!spacingMassRange[1].Contains(peakSpacings[1]))
            //                                {
            //                                    pair.peaks[channelIndex + 1, k] = null;
            //                                    peakCount--;

            //                                    // All inner spacings bad
            //                                    if (!spacingMassRange[2].Contains(peakSpacings[2]))
            //                                    {
            //                                        pair.peaks[channelIndex + 2, k] = null;
            //                                        peakCount--;
            //                                    }
            //                                }
            //                                // First & third inner spacings bad
            //                                if (!spacingMassRange[2].Contains(peakSpacings[2])) { }
            //                            }
            //                            // Second inner spacing bad
            //                            else if (!spacingMassRange[1].Contains(peakSpacings[1]))
            //                            {
            //                                // Second & third inner spacings bad
            //                                if (!spacingMassRange[2].Contains(peakSpacings[2]))
            //                                {
            //                                    pair.peaks[channelIndex + 2, k] = null;
            //                                    peakCount--;
            //                                }
            //                            }
            //                            // Third inner spacing bad
            //                            else if (!spacingMassRange[2].Contains(peakSpacings[2])) { }
            //                        }
            //                    }
            //                }
            //                // Incomplete peak sets
            //                else if (peakCount >= peaksNeeded && peakCount > 1)
            //                {
            //                    if (pair.peaks[0, k] != null)
            //                    {
            //                        // Peaks 1 & 2
            //                        if (pair.peaks[1, k] != null)
            //                        {
            //                            peakSpacings[0] = Mass.MassFromMz(pair.peaks[1, k].X - pair.peaks[0, k].X, charge);
            //                            if (!spacingMassRange[0].Contains(peakSpacings[0]))
            //                            {
            //                                pair.peaks[0, k] = null;
            //                                pair.peaks[1, k] = null;
            //                                peakCount -= 2;
            //                            }
            //                        }
            //                        // Peaks 1 & 3
            //                        else
            //                        {
            //                            peakSpacings[0] = Mass.MassFromMz(pair.peaks[2, k].X - pair.peaks[0, k].X, charge);
            //                            doubleSpacing = new MassRange(spacingMassRange[0].Mean + spacingMassRange[1].Mean, new MassTolerance(MassToleranceType.DA, 0.020 * numClusterLabels));
            //                            if (!doubleSpacing.Contains(peakSpacings[0]))
            //                            {
            //                                pair.peaks[0, k] = null;
            //                                pair.peaks[1, k] = null;
            //                                peakCount -= 2;
            //                            }
            //                        }
            //                    }
            //                    // Peaks 2 & 3
            //                    else
            //                    {
            //                        peakSpacings[0] = Mass.MassFromMz(pair.peaks[2, k].X - pair.peaks[1, k].X, charge);
            //                        if (!spacingMassRange[1].Contains(peakSpacings[0]))
            //                        {
            //                            pair.peaks[1, k] = null;
            //                            pair.peaks[2, k] = null;
            //                            peakCount -= 2;
            //                        }
            //                    }
            //                }

            //                if (peakCount < peaksNeeded)
            //                {
            //                    for (int j = channelIndex; j < channelIndex + numIsotopologues; j++)
            //                    {
            //                        pair.peaks[j, k] = null;
            //                    }
            //                }
            //            }
            //        }
            //    }
            //}

            //ILabeledPeak peak1;
            //ILabeledPeak peak2;
            //ILabeledPeak peak3;
            //ILabeledPeak peak4;
            //ILabeledPeak peak5;
            //ILabeledPeak peak6;

            //double spacing1;
            //double spacing2;
            //double spacing3;
            //double spacing4;
            //double spacing5;

            //double peak1NeutralMass;
            //double peak2NeutralMass;
            //double peak3NeutralMass;
            //double peak4NeutralMass;
            //double peak5NeutralMass;
            //double peak6NeutralMass;

            //if (numIsotopologues == 1)
            //{
            //    for (int i = 0; i < pairs.Count; i++)
            //    {
            //        pair = pairs[i];

            //        if (numChannels == 2)
            //        {
            //            for (int k = 0; k < numIsotopes; k++)
            //            {
            //                peakCount = 0;
            //                peak1 = pair.peaks[0, k];
            //                peak2 = pair.peaks[1, k];

            //                if (peak1 != null) peakCount++;
            //                if (peak2 != null) peakCount++;

            //                if (peakCount == 2)
            //                {
            //                    peak1NeutralMass = Mass.MassFromMz(peak1.X, charge);
            //                    peak2NeutralMass = Mass.MassFromMz(peak2.X, charge);

            //                    spacing1 = peak2NeutralMass - peak1NeutralMass;

            //                    //Remove pairs whose spacing is more than 10 mDa away from the calculated spacing
            //                    if (spacing1 < spacingMassRange[0].Minimum || spacing1 > spacingMassRange[0].Maximum)
            //                    {
            //                        //Console.WriteLine("bad spacing");
            //                        pair.peaks[0, k] = null;
            //                        peakCount--;
            //                        pair.peaks[1, k] = null;
            //                        peakCount--;
            //                    }
            //                }
            //            }
            //        }
            //        else if (numChannels == 3)
            //        {
            //            for (int k = 0; k < numIsotopes; k++)
            //            {
            //                peakCount = 0;
            //                peak1 = pair.peaks[0, k];
            //                peak2 = pair.peaks[1, k];
            //                peak3 = pair.peaks[2, k];

            //                if (peak1 != null) peakCount++;
            //                if (peak2 != null) peakCount++;
            //                if (peak3 != null) peakCount++;

            //                peakCount = 0;
            //                peak1 = pair.peaks[0, k];
            //                peak2 = pair.peaks[1, k];
            //                peak3 = pair.peaks[2, k];

            //                if (peak1 != null) peakCount++;
            //                if (peak2 != null) peakCount++;
            //                if (peak3 != null) peakCount++;

            //                if (peakCount == 3)
            //                {
            //                    peak1NeutralMass = Mass.MassFromMz(peak1.X, charge);
            //                    peak2NeutralMass = Mass.MassFromMz(peak2.X, charge);
            //                    peak3NeutralMass = Mass.MassFromMz(peak3.X, charge);

            //                    spacing1 = peak2NeutralMass - peak1NeutralMass;
            //                    spacing2 = peak3NeutralMass - peak2NeutralMass;

            //                    if (Form1.NOISEBANDCAP)
            //                    {
            //                        if (!spacingMassRange[0].Contains(spacing1))
            //                        {
            //                            pair.peaks[0, k] = null;
            //                            peakCount--;

            //                            if (!spacingMassRange[1].Contains(spacing2))
            //                            {
            //                                pair.peaks[1, k] = null;
            //                                peakCount--;
            //                                pair.peaks[2, k] = null;
            //                                peakCount--;
            //                            }
            //                        }
            //                        else if (!spacingMassRange[1].Contains(spacing2))
            //                        {
            //                            pair.peaks[2, k] = null;
            //                            peakCount--;
            //                        }

            //                        if (peakCount < peaksNeeded)
            //                        {
            //                            pair.peaks[0, k] = null;
            //                            pair.peaks[1, k] = null;
            //                            pair.peaks[2, k] = null;
            //                        }
            //                    }
            //                    else
            //                    {
            //                        if (!spacingMassRange[0].Contains(spacing1) || !spacingMassRange[1].Contains(spacing2))
            //                        {
            //                            pair.peaks[0, k] = null;
            //                            pair.peaks[1, k] = null;
            //                            pair.peaks[2, k] = null;
            //                        }
            //                    }
            //                }
            //                else if (Form1.NOISEBANDCAP && peakCount == 2)
            //                {
            //                    if (peak1 != null)
            //                    {
            //                        peak1NeutralMass = Mass.MassFromMz(peak1.X, charge);

            //                        // Peaks 1 & 2
            //                        if (peak2 != null)
            //                        {
            //                            peak2NeutralMass = Mass.MassFromMz(peak2.X, charge);
            //                            spacing1 = peak2NeutralMass - peak1NeutralMass;
            //                            if (!spacingMassRange[0].Contains(spacing1))
            //                            {
            //                                pair.peaks[0, k] = null;
            //                                pair.peaks[1, k] = null;
            //                            }
            //                        }
            //                        // Peaks 1 & 3
            //                        else if (peak3 != null)
            //                        {
            //                            peak2NeutralMass = Mass.MassFromMz(peak3.X, charge);
            //                            spacing1 = peak2NeutralMass - peak1NeutralMass;
            //                            double extendedTheoSpacing = spacingMassRange[0].Mean + spacingMassRange[1].Mean;
            //                            MassRange doubleSpacing = new MassRange(extendedTheoSpacing, new MassTolerance(MassToleranceType.DA, SPACINGPERCENTAGEERROR * extendedTheoSpacing));
            //                            if (!doubleSpacing.Contains(spacing1))
            //                            {
            //                                pair.peaks[0, k] = null;
            //                                pair.peaks[2, k] = null;
            //                            }
            //                        }
            //                    }
            //                    // Peaks 2 & 3
            //                    else if (peak2 != null)
            //                    {
            //                        peak1NeutralMass = Mass.MassFromMz(peak2.X, charge);
            //                        peak2NeutralMass = Mass.MassFromMz(peak3.X, charge);
            //                        spacing1 = peak2NeutralMass - peak1NeutralMass;
            //                        if (!spacingMassRange[1].Contains(spacing1))
            //                        {
            //                            pair.peaks[1, k] = null;
            //                            pair.peaks[2, k] = null;
            //                        }
            //                    }
            //                }
            //                else if (Form1.NOISEBANDCAP && peakCount < peaksNeeded)
            //                {
            //                    pair.peaks[0, k] = null;
            //                    pair.peaks[1, k] = null;
            //                    pair.peaks[2, k] = null;
            //                }

            //                if (peakCount < peaksNeeded)
            //                {
            //                    pair.peaks[0, k] = null;
            //                    pair.peaks[1, k] = null;
            //                    pair.peaks[2, k] = null;
            //                }
            //            }
            //        }
            //    }
            //}
            //if (numIsotopologues == 2)
            //{
            //    for (int i = 0; i < pairs.Count; i++)
            //    {
            //        pair = pairs[i];

            //        for (int j = 0; j < numChannels; j += 2)
            //        {
            //            for (int k = 0; k < numIsotopes; k++)
            //            {
            //                peakCount = 0;
            //                peak1 = pair.peaks[j, k];
            //                peak2 = pair.peaks[j + 1, k];

            //                if (peak1 != null) peakCount++;
            //                if (peak2 != null) peakCount++;

            //                // Only look at complete pairs
            //                if (peakCount == 2)
            //                {
            //                    peak1NeutralMass = Mass.MassFromMz(peak1.X, charge);
            //                    peak2NeutralMass = Mass.MassFromMz(peak2.X, charge);

            //                    spacing1 = peak2NeutralMass - peak1NeutralMass;

            //                    //Remove pairs whose spacing is more than 5 mDa away from the calculated spacing
            //                    if (spacing1 < spacingMassRange[0].Minimum || spacing1 > spacingMassRange[0].Maximum)
            //                    {
            //                        //Console.WriteLine("bad spacing");
            //                        pair.peaks[j, k] = null;
            //                        peakCount--;
            //                        pair.peaks[j + 1, k] = null;
            //                        peakCount--;
            //                    }
            //                }

            //                if (peakCount < peaksNeeded)
            //                {
            //                    pair.peaks[j, k] = null;
            //                    pair.peaks[j + 1, k] = null;
            //                }
            //            }
            //        }
            //    }
            //}
            //else if (numIsotopologues == 3)
            //{
            //    for (int i = 0; i < pairs.Count; i++)
            //    {
            //        pair = pairs[i];

            //        for (int j = 0; j < numChannels; j += 3)
            //        {
            //            for (int k = 0; k < numIsotopes; k++)
            //            {
            //                peakCount = 0;
            //                peak1 = pair.peaks[j, k];
            //                peak2 = pair.peaks[j + 1, k];
            //                peak3 = pair.peaks[j + 2, k];

            //                if (peak1 != null) peakCount++;
            //                if (peak2 != null) peakCount++;
            //                if (peak3 != null) peakCount++;

            //                // Only look at complete pairs
            //                if (peakCount == 3)
            //                {
            //                    peak1NeutralMass = Mass.MassFromMz(peak1.X, charge);
            //                    peak2NeutralMass = Mass.MassFromMz(peak2.X, charge);
            //                    peak3NeutralMass = Mass.MassFromMz(peak3.X, charge);

            //                    spacing1 = peak2NeutralMass - peak1NeutralMass;
            //                    spacing2 = peak3NeutralMass - peak2NeutralMass;

            //                    //Remove pairs whose spacing is more than 5 mDa away from the calculated spacing
            //                    if (Form1.NOISEBANDCAP)
            //                    {
            //                        if (!spacingMassRange[0].Contains(spacing1))
            //                        {
            //                            pair.peaks[j, k] = null;
            //                            peakCount--;

            //                            if (!spacingMassRange[1].Contains(spacing2))
            //                            {
            //                                pair.peaks[j + 1, k] = null;
            //                                peakCount--;
            //                                pair.peaks[j + 2, k] = null;
            //                                peakCount--;
            //                            }
            //                        }
            //                        else if (!spacingMassRange[1].Contains(spacing2))
            //                        {
            //                            pair.peaks[j + 2, k] = null;
            //                            peakCount--;
            //                        }

            //                        if (peakCount < peaksNeeded)
            //                        {
            //                            pair.peaks[j, k] = null;
            //                            pair.peaks[j + 1, k] = null;
            //                            pair.peaks[j + 2, k] = null;
            //                        }
            //                    }
            //                    else
            //                    {
            //                        if (!spacingMassRange[0].Contains(spacing1) || !spacingMassRange[1].Contains(spacing2))
            //                        {
            //                            pair.peaks[j, k] = null;
            //                            pair.peaks[j + 1, k] = null;
            //                            pair.peaks[j + 2, k] = null;
            //                        }
            //                    }
            //                }
            //                else if (Form1.NOISEBANDCAP && peakCount == 2)
            //                {
            //                    if (peak1 != null)
            //                    {
            //                        peak1NeutralMass = Mass.MassFromMz(peak1.X, charge);

            //                        // Peaks 1 & 2
            //                        if (peak2 != null)
            //                        {
            //                            peak2NeutralMass = Mass.MassFromMz(peak2.X, charge);
            //                            spacing1 = peak2NeutralMass - peak1NeutralMass;
            //                            if (!spacingMassRange[0].Contains(spacing1))
            //                            {
            //                                pair.peaks[j, k] = null;
            //                                pair.peaks[j + 1, k] = null;
            //                            }
            //                        }
            //                        // Peaks 1 & 3
            //                        else if (peak3 != null)
            //                        {
            //                            peak2NeutralMass = Mass.MassFromMz(peak3.X, charge);
            //                            spacing1 = peak2NeutralMass - peak1NeutralMass;
            //                            double extendedTheoSpacing = spacingMassRange[0].Mean + spacingMassRange[1].Mean;
            //                            MassRange doubleSpacing = new MassRange(extendedTheoSpacing, new MassTolerance(MassToleranceType.DA, SPACINGPERCENTAGEERROR * extendedTheoSpacing));
            //                            if (!doubleSpacing.Contains(spacing1))
            //                            {
            //                                pair.peaks[j, k] = null;
            //                                pair.peaks[j + 2, k] = null;
            //                            }
            //                        }
            //                    }
            //                    // Peaks 2 & 3
            //                    else if (peak2 != null)
            //                    {
            //                        peak1NeutralMass = Mass.MassFromMz(peak2.X, charge);
            //                        peak2NeutralMass = Mass.MassFromMz(peak3.X, charge);
            //                        spacing1 = peak2NeutralMass - peak1NeutralMass;
            //                        if (!spacingMassRange[1].Contains(spacing1))
            //                        {
            //                            pair.peaks[j + 1, k] = null;
            //                            pair.peaks[j + 2, k] = null;
            //                        }
            //                    }
            //                }
            //                else if (Form1.NOISEBANDCAP && peakCount < 1)
            //                {
            //                    pair.peaks[j, k] = null;
            //                    pair.peaks[j + 1, k] = null;
            //                    pair.peaks[j + 2, k] = null;
            //                }

            //                if (peakCount < peaksNeeded)
            //                {
            //                    pair.peaks[j, k] = null;
            //                    pair.peaks[j + 1, k] = null;
            //                    pair.peaks[j + 2, k] = null;
            //                }
            //            }
            //        }
            //    }
            //}
            //else if (numIsotopologues == 4)
            //{
            //    for (int i = 0; i < pairs.Count; i++)
            //    {
            //        pair = pairs[i];
            //        for (int j = 0; j < numChannels; j += 4)
            //        {
            //            for (int k = 0; k < numIsotopes; k++)
            //            {
            //                peakCount = 0;
            //                peak1 = pair.peaks[j, k];
            //                peak2 = pair.peaks[j + 1, k];
            //                peak3 = pair.peaks[j + 2, k];
            //                peak4 = pair.peaks[j + 3, k];

            //                if (peak1 != null) peakCount++;
            //                if (peak2 != null) peakCount++;
            //                if (peak3 != null) peakCount++;
            //                if (peak4 != null) peakCount++;

            //                // Only look at complete pairs
            //                if (peakCount == 4)
            //                {
            //                    peak1NeutralMass = Mass.MassFromMz(peak1.X, charge);
            //                    peak2NeutralMass = Mass.MassFromMz(peak2.X, charge);
            //                    peak3NeutralMass = Mass.MassFromMz(peak3.X, charge);
            //                    peak4NeutralMass = Mass.MassFromMz(peak4.X, charge);

            //                    spacing1 = peak2NeutralMass - peak1NeutralMass;
            //                    spacing2 = peak3NeutralMass - peak2NeutralMass;
            //                    spacing3 = peak4NeutralMass - peak3NeutralMass;

            //                    //Remove pairs whose spacing is more than 5 mDa away from the calculated spacing
            //                    if (Form1.NOISEBANDCAP)
            //                    {
            //                        // Spacing 1 is bad
            //                        if (!spacingMassRange[0].Contains(spacing1))
            //                        {
            //                            pair.peaks[j, k] = null;
            //                            peakCount--;

            //                            // Spacings 1 & 2 are bad
            //                            if (!spacingMassRange[1].Contains(spacing2))
            //                            {
            //                                pair.peaks[j + 1, k] = null;
            //                                peakCount--;

            //                                // Spacings 1, 2 & 3 are bad 
            //                                if (!spacingMassRange[2].Contains(spacing3))
            //                                {
            //                                    pair.peaks[j + 2, k] = null;
            //                                    peakCount--;
            //                                    pair.peaks[j + 3, k] = null;
            //                                    peakCount--;
            //                                }
            //                            }
            //                            // Spacings 1 & 3 are bad
            //                            else if (!spacingMassRange[2].Contains(spacing3))
            //                            {
            //                                pair.peaks[j + 3, k] = null;
            //                                peakCount--;
            //                            }
            //                        }
            //                        // Spacing 2 is bad
            //                        else if (!spacingMassRange[1].Contains(spacing2))
            //                        {
            //                            // Spacings 2 & 3 are bad
            //                            if (!spacingMassRange[2].Contains(spacing3))
            //                            {
            //                                pair.peaks[j + 2, k] = null;
            //                                peakCount--;
            //                                pair.peaks[j + 3, k] = null;
            //                                peakCount--;
            //                            }
            //                        }
            //                        // Spacing 3 is bad
            //                        else if (!spacingMassRange[2].Contains(spacing3))
            //                        {
            //                            pair.peaks[j + 3, k] = null;
            //                            peakCount--;
            //                        }
            //                    }
            //                    else
            //                    {
            //                        if (!spacingMassRange[0].Contains(spacing1) || !spacingMassRange[1].Contains(spacing2) || !spacingMassRange[2].Contains(spacing3))
            //                        {
            //                            pair.peaks[j, k] = null;
            //                            pair.peaks[j + 1, k] = null;
            //                            pair.peaks[j + 2, k] = null;
            //                            pair.peaks[j + 3, k] = null;
            //                        }
            //                    }
            //                }
            //                else if (Form1.NOISEBANDCAP && peakCount == 3)
            //                {
            //                    if (peak1 != null)
            //                    {
            //                        peak1NeutralMass = Mass.MassFromMz(peak1.X, charge);

            //                        if (peak2 != null)
            //                        {
            //                            peak2NeutralMass = Mass.MassFromMz(peak2.X, charge);
            //                            spacing1 = peak2NeutralMass - peak1NeutralMass;

            //                            // Peaks 1, 2 & 3
            //                            if (peak3 != null)
            //                            {
            //                                peak3NeutralMass = Mass.MassFromMz(peak3.X, charge);
            //                                spacing2 = peak3NeutralMass - peak2NeutralMass;
            //                                spacing3 = peak3NeutralMass - peak1NeutralMass;
            //                                double outerTheoSpacing = spacingMassRange[0].Mean + spacingMassRange[1].Mean;
            //                                MassRange outerSpacing = new MassRange(outerTheoSpacing, new MassTolerance(MassToleranceType.DA, SPACINGPERCENTAGEERROR * outerTheoSpacing));

            //                                // Outer spacing is bad
            //                                if (!outerSpacing.Contains(outerTheoSpacing))
            //                                {
            //                                    // Both inner spacings are bad
            //                                    if (!spacingMassRange[0].Contains(spacing1) && !spacingMassRange[1].Contains(spacing2))
            //                                    {
            //                                        pair.peaks[j, k] = null;
            //                                        peakCount--;
            //                                        pair.peaks[j + 1, k] = null;
            //                                        peakCount--;
            //                                        pair.peaks[j + 2, k] = null;
            //                                        peakCount--;
            //                                    }
            //                                    // First inner spacing is bad
            //                                    else if (!spacingMassRange[0].Contains(spacing1))
            //                                    {
            //                                        pair.peaks[j, k] = null;
            //                                        peakCount--;
            //                                    }
            //                                    // Second inner spacing is bad
            //                                    else if (!spacingMassRange[1].Contains(spacing2))
            //                                    {
            //                                        pair.peaks[j + 2, k] = null;
            //                                        peakCount--;
            //                                    }
            //                                }
            //                                // Outer spacing is good
            //                                else
            //                                {
            //                                    // One inner spacing is bad
            //                                    if (!spacingMassRange[0].Contains(spacing1) || !spacingMassRange[1].Contains(spacing2))
            //                                    {
            //                                        pair.peaks[j + 1, k] = null;
            //                                        peakCount--;
            //                                    }
            //                                }
            //                            }
            //                            // Peaks 1, 2 & 4
            //                            else if (peak4 != null)
            //                            {
            //                                peak4NeutralMass = Mass.MassFromMz(peak4.X, charge);
            //                                spacing2 = peak4NeutralMass - peak2NeutralMass;
            //                                spacing3 = peak4NeutralMass - peak1NeutralMass;
            //                                double extendedTheoSpacing = spacingMassRange[1].Mean + spacingMassRange[2].Mean;
            //                                MassRange doubleSpacing = new MassRange(extendedTheoSpacing, new MassTolerance(MassToleranceType.DA, SPACINGPERCENTAGEERROR * extendedTheoSpacing));
            //                                double outerTheoSpacing = spacingMassRange[0].Mean + spacingMassRange[1].Mean + spacingMassRange[2].Mean;
            //                                MassRange outerSpacing = new MassRange(outerTheoSpacing, new MassTolerance(MassToleranceType.DA, SPACINGPERCENTAGEERROR * outerTheoSpacing));

            //                                // Outer spacing is bad
            //                                if (!outerSpacing.Contains(outerTheoSpacing))
            //                                {
            //                                    // Both inner spacings are bad
            //                                    if (!spacingMassRange[0].Contains(spacing1) && !doubleSpacing.Contains(spacing2))
            //                                    {
            //                                        pair.peaks[j, k] = null;
            //                                        peakCount--;
            //                                        pair.peaks[j + 1, k] = null;
            //                                        peakCount--;
            //                                        pair.peaks[j + 3, k] = null;
            //                                        peakCount--;
            //                                    }
            //                                    // First inner spacing is bad
            //                                    else if (!spacingMassRange[0].Contains(spacing1))
            //                                    {
            //                                        pair.peaks[j, k] = null;
            //                                        peakCount--;
            //                                    }
            //                                    // Second inner spacing is bad
            //                                    else if (!doubleSpacing.Contains(spacing2))
            //                                    {
            //                                        pair.peaks[j + 3, k] = null;
            //                                        peakCount--;
            //                                    }
            //                                }
            //                                // Outer spacing is good
            //                                else
            //                                {
            //                                    // One inner spacing is bad
            //                                    if (!spacingMassRange[0].Contains(spacing1) || !doubleSpacing.Contains(spacing2))
            //                                    {
            //                                        pair.peaks[j + 1, k] = null;
            //                                        peakCount--;
            //                                    }
            //                                }
            //                            }
            //                        }
            //                        else if (peak3 != null)
            //                        {
            //                            // Peaks 1, 3 & 4
            //                            peak3NeutralMass = Mass.MassFromMz(peak3.X, charge);
            //                            peak4NeutralMass = Mass.MassFromMz(peak4.X, charge);

            //                            spacing1 = peak3NeutralMass - peak1NeutralMass;
            //                            spacing2 = peak4NeutralMass - peak3NeutralMass;
            //                            spacing3 = peak4NeutralMass - peak1NeutralMass;

            //                            double extendedTheoSpacing = spacingMassRange[0].Mean + spacingMassRange[1].Mean;
            //                            MassRange doubleSpacing = new MassRange(extendedTheoSpacing, new MassTolerance(MassToleranceType.DA, SPACINGPERCENTAGEERROR * extendedTheoSpacing));
            //                            double outerTheoSpacing = spacingMassRange[0].Mean + spacingMassRange[1].Mean + spacingMassRange[2].Mean;
            //                            MassRange outerSpacing = new MassRange(outerTheoSpacing, new MassTolerance(MassToleranceType.DA, SPACINGPERCENTAGEERROR * outerTheoSpacing));

            //                            // Outer spacing is bad
            //                            if (!outerSpacing.Contains(spacing3))
            //                            {
            //                                // Both inner spacings are bad
            //                                if (!doubleSpacing.Contains(spacing1) && spacingMassRange[2].Contains(spacing2))
            //                                {
            //                                    pair.peaks[j, k] = null;
            //                                    peakCount--;
            //                                    pair.peaks[j + 2, k] = null;
            //                                    peakCount--;
            //                                    pair.peaks[j + 3, k] = null;
            //                                    peakCount--;
            //                                }

            //                                // First inner spacing is bad
            //                                else if (!doubleSpacing.Contains(spacing1))
            //                                {
            //                                    pair.peaks[j, k] = null;
            //                                    peakCount--;
            //                                }

            //                                // Second inner spacing is bad
            //                                else if (spacingMassRange[2].Contains(spacing2))
            //                                {
            //                                    pair.peaks[j + 3, k] = null;
            //                                    peakCount--;
            //                                }
            //                            }
            //                            // Outer spacing is good
            //                            else
            //                            {
            //                                // One inner spacing is bad
            //                                if (!doubleSpacing.Contains(spacing1) || spacingMassRange[2].Contains(spacing2))
            //                                {
            //                                    pair.peaks[j + 2, k] = null;
            //                                    peakCount--;
            //                                }
            //                            }
            //                        }
            //                    }
            //                    // Peaks 2, 3 & 4
            //                    else if (peak2 != null)
            //                    {
            //                        peak2NeutralMass = Mass.MassFromMz(peak2.X, charge);
            //                        peak3NeutralMass = Mass.MassFromMz(peak3.X, charge);
            //                        peak4NeutralMass = Mass.MassFromMz(peak4.X, charge);

            //                        spacing1 = peak3NeutralMass - peak2NeutralMass;
            //                        spacing2 = peak4NeutralMass - peak3NeutralMass;
            //                        spacing3 = peak4NeutralMass - peak2NeutralMass;

            //                        double outerTheoSpacing = spacingMassRange[1].Mean + spacingMassRange[2].Mean;
            //                        MassRange outerSpacing = new MassRange(outerTheoSpacing, new MassTolerance(MassToleranceType.DA, SPACINGPERCENTAGEERROR * outerTheoSpacing));

            //                        // Outer spacing is bad
            //                        if (!outerSpacing.Contains(spacing3))
            //                        {
            //                            // Both inner spacings are bad
            //                            if (!spacingMassRange[1].Contains(spacing1) && !spacingMassRange[2].Contains(spacing2))
            //                            {
            //                                pair.peaks[j + 1, k] = null;
            //                                peakCount--;
            //                                pair.peaks[j + 2, k] = null;
            //                                peakCount--;
            //                                pair.peaks[j + 3, k] = null;
            //                                peakCount--;
            //                            }

            //                            // First inner spacing is bad
            //                            else if (!spacingMassRange[1].Contains(spacing1))
            //                            {
            //                                pair.peaks[j + 1, k] = null;
            //                                peakCount--;
            //                            }

            //                            // Second inner spacing is bad
            //                            else if (!spacingMassRange[2].Contains(spacing2))
            //                            {
            //                                pair.peaks[j + 3, k] = null;
            //                                peakCount--;
            //                            }
            //                        }
            //                        // Outer spacing is good
            //                        else
            //                        {
            //                            // One inner spacing is bad
            //                            if (!spacingMassRange[1].Contains(spacing1) || !spacingMassRange[2].Contains(spacing2))
            //                            {
            //                                pair.peaks[j + 2, k] = null;
            //                                peakCount--;
            //                            }
            //                        }
            //                    }
            //                }
            //                else if (Form1.NOISEBANDCAP && peakCount == 2)
            //                {
            //                    if (peak1 != null)
            //                    {
            //                        peak1NeutralMass = Mass.MassFromMz(peak1.X, charge);

            //                        // Peaks 1 & 2
            //                        if (peak2 != null)
            //                        {
            //                            peak2NeutralMass = Mass.MassFromMz(peak2.X, charge);
            //                            spacing1 = peak2NeutralMass - peak1NeutralMass;
            //                            if (!spacingMassRange[0].Contains(spacing1))
            //                            {
            //                                pair.peaks[j, k] = null;
            //                                peakCount--;
            //                                pair.peaks[j + 1, k] = null;
            //                                peakCount--;
            //                            }
            //                        }
            //                        // Peaks 1 & 3
            //                        else if (peak3 != null)
            //                        {
            //                            peak2NeutralMass = Mass.MassFromMz(peak3.X, charge);
            //                            spacing1 = peak2NeutralMass - peak1NeutralMass;
            //                            double extendedTheoMass = spacingMassRange[0].Mean + spacingMassRange[1].Mean;
            //                            MassRange doubleSpacing = new MassRange(extendedTheoMass, new MassTolerance(MassToleranceType.DA, SPACINGPERCENTAGEERROR * extendedTheoMass));
            //                            if (!doubleSpacing.Contains(spacing1))
            //                            {
            //                                pair.peaks[j, k] = null;
            //                                peakCount--;
            //                                pair.peaks[j + 2, k] = null;
            //                                peakCount--;
            //                            }
            //                        }
            //                        // Peaks 1 & 4
            //                        else if (peak4 != null)
            //                        {
            //                            peak2NeutralMass = Mass.MassFromMz(peak4.X, charge);
            //                            spacing1 = peak2NeutralMass - peak1NeutralMass;
            //                            double extendedTheoMass = spacingMassRange[0].Mean + spacingMassRange[1].Mean + spacingMassRange[2].Mean;
            //                            MassRange doubleSpacing = new MassRange(extendedTheoMass, new MassTolerance(MassToleranceType.DA, SPACINGPERCENTAGEERROR * extendedTheoMass));
            //                            if (!doubleSpacing.Contains(spacing1))
            //                            {
            //                                pair.peaks[j, k] = null;
            //                                peakCount--;
            //                                pair.peaks[j + 3, k] = null;
            //                                peakCount--;
            //                            }
            //                        }
            //                    }
            //                    else if (peak2 != null)
            //                    {
            //                        peak1NeutralMass = Mass.MassFromMz(peak2.X, charge);

            //                        // Peaks 2 & 3
            //                        if (peak3 != null)
            //                        {
            //                            peak2NeutralMass = Mass.MassFromMz(peak3.X, charge);
            //                            spacing1 = peak2NeutralMass - peak1NeutralMass;
            //                            if (!spacingMassRange[1].Contains(spacing1))
            //                            {
            //                                pair.peaks[j + 1, k] = null;
            //                                peakCount--;
            //                                pair.peaks[j + 2, k] = null;
            //                                peakCount--;
            //                            }
            //                        }
            //                        // Peaks 2 & 4
            //                        else if (peak4 != null)
            //                        {
            //                            peak2NeutralMass = Mass.MassFromMz(peak4.X, charge);
            //                            spacing1 = peak2NeutralMass - peak1NeutralMass;
            //                            double extendedTheoSpacing = spacingMassRange[1].Mean + spacingMassRange[2].Mean;
            //                            MassRange doubleSpacing = new MassRange(extendedTheoSpacing, new MassTolerance(MassToleranceType.DA, SPACINGPERCENTAGEERROR * extendedTheoSpacing));
            //                            if (!doubleSpacing.Contains(spacing1))
            //                            {
            //                                pair.peaks[j + 1, k] = null;
            //                                peakCount--;
            //                                pair.peaks[j + 3, k] = null;
            //                                peakCount--;
            //                            }
            //                        }
            //                    }
            //                    // Peaks 3 & 4
            //                    else if (peak3 != null)
            //                    {
            //                        peak1NeutralMass = Mass.MassFromMz(peak3.X, charge);
            //                        peak2NeutralMass = Mass.MassFromMz(peak4.X, charge);
            //                        spacing1 = peak2NeutralMass - peak1NeutralMass;
            //                        if (!spacingMassRange[2].Contains(spacing1))
            //                        {
            //                            pair.peaks[j + 2, k] = null;
            //                            peakCount--;
            //                            pair.peaks[j + 3, k] = null;
            //                            peakCount--;
            //                        }
            //                    }
            //                }
            //                else if (Form1.NOISEBANDCAP && peakCount < 2)
            //                {
            //                    pair.peaks[j, k] = null;
            //                    pair.peaks[j + 1, k] = null;
            //                    pair.peaks[j + 2, k] = null;
            //                    pair.peaks[j + 3, k] = null;
            //                }

            //                if (peakCount < peaksNeeded)
            //                {
            //                    pair.peaks[j, k] = null;
            //                    pair.peaks[j + 1, k] = null;
            //                    pair.peaks[j + 2, k] = null;
            //                    pair.peaks[j + 3, k] = null;
            //                }
            //            }
            //        }
            //    }
            //}
            //else if (numIsotopologues == 6)
            //{
            //    for (int i = 0; i < pairs.Count; i++)
            //    {
            //        pair = pairs[i];
            //        for (int j = 0; j < numChannels; j += 6)
            //        {
            //            for (int k = 0; k < numIsotopes; k++)
            //            {
            //                peakCount = 0;
            //                peak1 = pair.peaks[j, k];
            //                peak2 = pair.peaks[j + 1, k];
            //                peak3 = pair.peaks[j + 2, k];
            //                peak4 = pair.peaks[j + 3, k];
            //                peak5 = pair.peaks[j + 4, k];
            //                peak6 = pair.peaks[j + 5, k];

            //                if (peak1 != null) peakCount++;
            //                if (peak2 != null) peakCount++;
            //                if (peak3 != null) peakCount++;
            //                if (peak4 != null) peakCount++;
            //                if (peak5 != null) peakCount++;
            //                if (peak6 != null) peakCount++;

            //                // Only look at complete pairs
            //                if (peakCount == 6)
            //                {
            //                    peak1NeutralMass = Mass.MassFromMz(peak1.X, charge);
            //                    peak2NeutralMass = Mass.MassFromMz(peak2.X, charge);
            //                    peak3NeutralMass = Mass.MassFromMz(peak3.X, charge);
            //                    peak4NeutralMass = Mass.MassFromMz(peak4.X, charge);
            //                    peak5NeutralMass = Mass.MassFromMz(peak5.X, charge);
            //                    peak6NeutralMass = Mass.MassFromMz(peak6.X, charge);

            //                    spacing1 = peak2NeutralMass - peak1NeutralMass;
            //                    spacing2 = peak3NeutralMass - peak2NeutralMass;
            //                    spacing3 = peak4NeutralMass - peak3NeutralMass;
            //                    spacing4 = peak5NeutralMass - peak4NeutralMass;
            //                    spacing5 = peak6NeutralMass - peak5NeutralMass;

            //                    //Remove pairs whose spacing is more than 5 mDa away from the calculated spacing
            //                    if (Form1.NOISEBANDCAP)
            //                    {
            //                        // Spacing 1 is bad
            //                        if (!spacingMassRange[0].Contains(spacing1))
            //                        {
            //                            pair.peaks[j, k] = null;
            //                            peakCount--;

            //                            // Spacings 1 & 2 are bad
            //                            if (!spacingMassRange[1].Contains(spacing2))
            //                            {
            //                                pair.peaks[j + 1, k] = null;
            //                                peakCount--;

            //                                // Spacings 1,2 & 3 are bad
            //                                if (!spacingMassRange[2].Contains(spacing3))
            //                                {
            //                                    pair.peaks[j + 2, k] = null;
            //                                    peakCount--;

            //                                    // Spacings 1,2,3 & 4 are bad
            //                                    if (!spacingMassRange[3].Contains(spacing4))
            //                                    {
            //                                        pair.peaks[j + 3, k] = null;
            //                                        peakCount--;

            //                                        // Spacings 1,2,3,4 & 5 are bad
            //                                        if (!spacingMassRange[4].Contains(spacing5))
            //                                        {
            //                                            pair.peaks[j + 4, k] = null;
            //                                            peakCount--;
            //                                            pair.peaks[j + 5, k] = null;
            //                                            peakCount--;
            //                                        }
            //                                    }
            //                                    // Spacings 1,2,3 & 5 are bad
            //                                    else if (!spacingMassRange[4].Contains(spacing5))
            //                                    {
            //                                        pair.peaks[j + 5, k] = null;
            //                                        peakCount--;
            //                                    }
            //                                }
            //                                // Spacings 1,2 & 4 are bad
            //                                else if (!spacingMassRange[3].Contains(spacing4))
            //                                {
            //                                    // Spacings 1,2,4 & 5 are bad
            //                                    if (!spacingMassRange[4].Contains(spacing5))
            //                                    {
            //                                        pair.peaks[j + 4, k] = null;
            //                                        peakCount--;
            //                                        pair.peaks[j + 5, k] = null;
            //                                        peakCount--;
            //                                    }
            //                                }
            //                            }
            //                            // Spacings 1 & 3 are bad
            //                            else if (!spacingMassRange[2].Contains(spacing3))
            //                            {
            //                                // Spacings 1,3 & 4 are bad
            //                                if (!spacingMassRange[3].Contains(spacing4))
            //                                {
            //                                    pair.peaks[j + 3, k] = null;
            //                                    peakCount--;

            //                                    // Spacings 1,3,4 & 5 are bad
            //                                    if (!spacingMassRange[4].Contains(spacing5))
            //                                    {
            //                                        pair.peaks[j + 4, k] = null;
            //                                        peakCount--;
            //                                        pair.peaks[j + 5, k] = null;
            //                                        peakCount--;
            //                                    }
            //                                }
            //                                // Spacings 1,3 & 5 are bad
            //                                else if (!spacingMassRange[4].Contains(spacing5))
            //                                {
            //                                    pair.peaks[j + 5, k] = null;
            //                                    peakCount--;
            //                                }
            //                            }
            //                            // Spacings 1 & 4 are bad
            //                            else if (!spacingMassRange[3].Contains(spacing4))
            //                            {
            //                                // Spacings 1,4 & 5 are bad
            //                                if (!spacingMassRange[4].Contains(spacing5))
            //                                {
            //                                    pair.peaks[j + 4, k] = null;
            //                                    peakCount--;
            //                                    pair.peaks[j + 5, k] = null;
            //                                    peakCount--;
            //                                }
            //                            }
            //                            // Spacings 1 & 5 are bad
            //                            else if (!spacingMassRange[4].Contains(spacing5))
            //                            {
            //                                pair.peaks[j + 5, k] = null;
            //                                peakCount--;
            //                            }
            //                        }
            //                        // Spacing 2 is bad
            //                        else if (!spacingMassRange[1].Contains(spacing2))
            //                        {
            //                            // Spacings 2 & 3 are bad
            //                            if (!spacingMassRange[2].Contains(spacing3))
            //                            {
            //                                pair.peaks[j + 2, k] = null;
            //                                peakCount--;
            //                                // Spacings 2,3 & 4 are bad
            //                                if (!spacingMassRange[3].Contains(spacing4))
            //                                {
            //                                    pair.peaks[j + 3, k] = null;
            //                                    peakCount--;
            //                                    // Spacings 2,3,4 & 5 are bad
            //                                    if (!spacingMassRange[4].Contains(spacing5))
            //                                    {
            //                                        pair.peaks[j + 4, k] = null;
            //                                        peakCount--;
            //                                        pair.peaks[j + 5, k] = null;
            //                                        peakCount--;
            //                                    }
            //                                }
            //                                // Spacings 2,3 & 5 are bad
            //                                else if (!spacingMassRange[4].Contains(spacing5))
            //                                {
            //                                    pair.peaks[j + 5, k] = null;
            //                                    peakCount--;
            //                                }
            //                            }
            //                            // Spacings 2 & 4 are bad
            //                            else if (!spacingMassRange[3].Contains(spacing4))
            //                            {
            //                                // Spacings 2,4 & 5 are bad
            //                                if (!spacingMassRange[4].Contains(spacing5))
            //                                {
            //                                    pair.peaks[j + 4, k] = null;
            //                                    peakCount--;
            //                                    pair.peaks[j + 5, k] = null;
            //                                    peakCount--;
            //                                }
            //                            }
            //                            // Spacings 2 & 5 are bad
            //                            else if (!spacingMassRange[4].Contains(spacing5))
            //                            {
            //                                pair.peaks[j + 5, k] = null;
            //                                peakCount--;
            //                            }
            //                        }
            //                        // Spacing 3 is bad
            //                        else if (!spacingMassRange[2].Contains(spacing3))
            //                        {
            //                            // Spacings 3 & 4 are bad
            //                            if (!spacingMassRange[3].Contains(spacing4))
            //                            {
            //                                pair.peaks[j + 3, k] = null;
            //                                peakCount--;
            //                                // Spacings 3,4 & 5 are bad
            //                                if (!spacingMassRange[4].Contains(spacing5))
            //                                {
            //                                    pair.peaks[j + 4, k] = null;
            //                                    peakCount--;
            //                                    pair.peaks[j + 5, k] = null;
            //                                    peakCount--;
            //                                }
            //                            }
            //                            // Spacings 3 & 5 are bad
            //                            else if (!spacingMassRange[4].Contains(spacing5))
            //                            {
            //                                pair.peaks[j + 5, k] = null;
            //                                peakCount--;
            //                            }
            //                        }
            //                        // Spacing 4 is bad
            //                        else if (!spacingMassRange[3].Contains(spacing4))
            //                        {
            //                            // Spacings 4 & 5 are bad
            //                            if (!spacingMassRange[4].Contains(spacing5))
            //                            {
            //                                pair.peaks[j + 4, k] = null;
            //                                peakCount--;
            //                                pair.peaks[j + 5, k] = null;
            //                                peakCount--;
            //                            }
            //                        }
            //                        // Spacing 5 is bad
            //                        else if (!spacingMassRange[4].Contains(spacing5))
            //                        {
            //                            pair.peaks[j + 5, k] = null;
            //                            peakCount--;
            //                        }
            //                    }
            //                    else
            //                    {
            //                        if (!spacingMassRange[0].Contains(spacing1) || !spacingMassRange[1].Contains(spacing2) || !spacingMassRange[2].Contains(spacing3) || !spacingMassRange[3].Contains(spacing4) || !spacingMassRange[4].Contains(spacing5))
            //                        {
            //                            pair.peaks[j + 1, k] = null;
            //                            pair.peaks[j + 2, k] = null;
            //                            pair.peaks[j + 3, k] = null;
            //                            pair.peaks[j + 4, k] = null;
            //                            pair.peaks[j + 5, k] = null;
            //                        }
            //                    }
            //                }
            //                else if (Form1.NOISEBANDCAP && peakCount == 5)
            //                {
            //                    if (peak1 != null)
            //                    {
            //                        peak1NeutralMass = Mass.MassFromMz(peak1.X, charge);

            //                        if (peak2 != null)
            //                        {
            //                            peak2NeutralMass = Mass.MassFromMz(peak2.X, charge);
            //                            spacing1 = peak2NeutralMass - peak1NeutralMass;

            //                            if (peak3 != null)
            //                            {
            //                                peak3NeutralMass = Mass.MassFromMz(peak3.X, charge);
            //                                spacing2 = peak3NeutralMass - peak2NeutralMass;

            //                                if (peak4 != null)
            //                                {
            //                                    peak4NeutralMass = Mass.MassFromMz(peak4.X, charge);
            //                                    spacing3 = peak4NeutralMass - peak3NeutralMass;

            //                                    // Peaks 1,2,3,4 & 5
            //                                    if (peak5 != null)
            //                                    {
            //                                        peak5NeutralMass = Mass.MassFromMz(peak5.X, charge);
            //                                        spacing4 = peak5NeutralMass - peak4NeutralMass;
            //                                    }
            //                                    // Peaks 1,2,3,4 & 6
            //                                    else if (peak6 != null)
            //                                    {
            //                                        peak5NeutralMass = Mass.MassFromMz(peak6.X, charge);
            //                                        spacing4 = peak5NeutralMass - peak4NeutralMass;
            //                                        double extendedTheoSpacing = spacingMassRange[3].Mean + spacingMassRange[4].Mean;
            //                                        MassRange doubleSpacing = new MassRange(extendedTheoSpacing, new MassTolerance(MassToleranceType.DA, SPACINGPERCENTAGEERROR * extendedTheoSpacing));
            //                                    }
            //                                }
            //                            }
            //                            // Peaks 1,2,4,5 & 6
            //                            else if (peak4 != null)
            //                            {
            //                                peak3NeutralMass = Mass.MassFromMz(peak4.X, charge);
            //                                spacing2 = peak3NeutralMass - peak2NeutralMass;
            //                                double extendedTheoSpacing = spacingMassRange[1].Mean + spacingMassRange[2].Mean;
            //                                MassRange doubleSpacing = new MassRange(extendedTheoSpacing, new MassTolerance(MassToleranceType.DA, SPACINGPERCENTAGEERROR * extendedTheoSpacing));
            //                            }
            //                        }
            //                        // Peaks 1,3,4,5 & 6
            //                        else if (peak3 != null)
            //                        {
            //                            peak2NeutralMass = Mass.MassFromMz(peak3.X, charge);
            //                            spacing1 = peak2NeutralMass - peak1NeutralMass;
            //                            double extendedTheoSpacing = spacingMassRange[0].Mean + spacingMassRange[1].Mean;
            //                            MassRange doubleSpacing = new MassRange(extendedTheoSpacing, new MassTolerance(MassToleranceType.DA, SPACINGPERCENTAGEERROR * extendedTheoSpacing));
            //                        }
            //                    }
            //                    // Peaks 2,3,4,5 & 6
            //                    else if (peak2 != null)
            //                    {
            //                        peak1NeutralMass = Mass.MassFromMz(peak2.X, charge);
            //                        peak2NeutralMass = Mass.MassFromMz(peak3.X, charge);
            //                        peak3NeutralMass = Mass.MassFromMz(peak4.X, charge);
            //                        peak4NeutralMass = Mass.MassFromMz(peak5.X, charge);
            //                        peak5NeutralMass = Mass.MassFromMz(peak6.X, charge);

            //                        spacing1 = peak2NeutralMass - peak1NeutralMass;
            //                        spacing2 = peak3NeutralMass - peak2NeutralMass;
            //                        spacing3 = peak4NeutralMass - peak3NeutralMass;
            //                        spacing4 = peak5NeutralMass - peak4NeutralMass;

            //                        // Spacing 1 is bad
            //                        if (!spacingMassRange[1].Contains(spacing1))
            //                        {
            //                            pair.peaks[j + 1, k] = null;
            //                            peakCount--;

            //                            // Spacings 1 & 2 are bad
            //                            if (!spacingMassRange[2].Contains(spacing2))
            //                            {
            //                                pair.peaks[j + 2, k] = null;
            //                                peakCount--;

            //                                // Spacings 1,2 & 3 are bad
            //                                if (!spacingMassRange[3].Contains(spacing3))
            //                                {
            //                                    pair.peaks[j + 3, k] = null;
            //                                    peakCount--;

            //                                    // Spacings 1,2,3 & 4 are bad
            //                                    if (!spacingMassRange[4].Contains(spacing4))
            //                                    {
            //                                        pair.peaks[j + 4, k] = null;
            //                                        peakCount--;
            //                                        pair.peaks[j + 5, k] = null;
            //                                        peakCount--;
            //                                    }
            //                                }
            //                                // Spacings 1,2 & 4 are bad
            //                                else if (!spacingMassRange[4].Contains(spacing4))
            //                                {
            //                                    pair.peaks[j + 5, k] = null;
            //                                    peakCount--;
            //                                }
            //                            }
            //                            // Spacings 1 & 3 are bad
            //                            else if (!spacingMassRange[3].Contains(spacing3))
            //                            {
            //                                // Spacings 1,3 & 4 are bad
            //                                if (!spacingMassRange[4].Contains(spacing4))
            //                                {
            //                                    pair.peaks[j + 4, k] = null;
            //                                    peakCount--;
            //                                    pair.peaks[j + 5, k] = null;
            //                                    peakCount--;
            //                                }
            //                            }

            //                            // Spacings 1 & 4 are bad
            //                            else if (!spacingMassRange[4].Contains(spacing4))
            //                            {
            //                                pair.peaks[j + 5, k] = null;
            //                                peakCount--;
            //                            }
            //                        }
            //                        // Spacing 2 is bad
            //                        else if (!spacingMassRange[2].Contains(spacing2))
            //                        {
            //                            // Spacings 2 & 3 are bad
            //                            if (!spacingMassRange[3].Contains(spacing3))
            //                            {
            //                                pair.peaks[j + 3, k] = null;
            //                                peakCount--;
            //                                // Spacings 2,3 & 4 are bad
            //                                if (!spacingMassRange[4].Contains(spacing4))
            //                                {
            //                                    pair.peaks[j + 4, k] = null;
            //                                    peakCount--;
            //                                    pair.peaks[j + 5, k] = null;
            //                                    peakCount--;
            //                                }
            //                            }
            //                            // Spacings 2 & 4 are bad
            //                            else if (!spacingMassRange[4].Contains(spacing4))
            //                            {
            //                                pair.peaks[j + 5, k] = null;
            //                                peakCount--;
            //                            }
            //                        }
            //                        // Spacing 3 is bad
            //                        else if (!spacingMassRange[3].Contains(spacing3))
            //                        {
            //                            // Spacings 3 & 4 are bad
            //                            if (!spacingMassRange[4].Contains(spacing4))
            //                            {
            //                                pair.peaks[j + 4, k] = null;
            //                                peakCount--;
            //                                pair.peaks[j + 5, k] = null;
            //                                peakCount--;
            //                            }
            //                        }
            //                        // Spacing 4 is bad
            //                        else if (!spacingMassRange[4].Contains(spacing4))
            //                        {
            //                            pair.peaks[j + 5, k] = null;
            //                            peakCount--;
            //                        }
            //                    }
            //                }
            //                else if (Form1.NOISEBANDCAP && peakCount == 4)
            //                {
            //                    if (peak1 != null)
            //                    {
            //                        peak1NeutralMass = Mass.MassFromMz(peak1.X, charge);

            //                        if (peak2 != null)
            //                        {
            //                            peak2NeutralMass = Mass.MassFromMz(peak2.X, charge);
            //                            spacing1 = peak2NeutralMass - peak1NeutralMass;

            //                            if (peak3 != null)
            //                            {
            //                                peak3NeutralMass = Mass.MassFromMz(peak3.X, charge);
            //                                spacing2 = peak3NeutralMass - peak2NeutralMass;

            //                                // Peaks 1,2,3 & 4
            //                                if (peak4 != null)
            //                                {
            //                                    peak4NeutralMass = Mass.MassFromMz(peak4.X, charge);
            //                                    spacing3 = peak4NeutralMass - peak3NeutralMass;

            //                                    // Spacing 1 is bad
            //                                    if (!spacingMassRange[0].Contains(spacing1))
            //                                    {
            //                                        pair.peaks[j, k] = null;
            //                                        peakCount--;

            //                                        // Spacings 1 & 2 are bad
            //                                        if (!spacingMassRange[1].Contains(spacing2))
            //                                        {
            //                                            pair.peaks[j + 1, k] = null;
            //                                            peakCount--;

            //                                            // Spacings 1,2 & 3 are bad
            //                                            if (!spacingMassRange[2].Contains(spacing3))
            //                                            {
            //                                                pair.peaks[j + 2, k] = null;
            //                                                peakCount--;
            //                                                pair.peaks[j + 3, k] = null;
            //                                                peakCount--;
            //                                            }
            //                                        }
            //                                        // Spacings 1 & 3 are bad
            //                                        else if (!spacingMassRange[2].Contains(spacing3))
            //                                        {
            //                                            pair.peaks[j + 3, k] = null;
            //                                            peakCount--;
            //                                        }
            //                                    }
            //                                    // Spacing 2 is bad
            //                                    else if (!spacingMassRange[1].Contains(spacing2))
            //                                    {
            //                                        // Spacings 2 & 3 are bad
            //                                        if (!spacingMassRange[2].Contains(spacing3))
            //                                        {
            //                                            pair.peaks[j + 2, k] = null;
            //                                            peakCount--;
            //                                            pair.peaks[j + 3, k] = null;
            //                                            peakCount--;
            //                                        }
            //                                    }
            //                                    // Spacing 3 is bad
            //                                    else if (!spacingMassRange[2].Contains(spacing3))
            //                                    {
            //                                        pair.peaks[j + 3, k] = null;
            //                                        peakCount--;
            //                                    }
            //                                }
            //                                // Peaks 1,2,3 & 5
            //                                else if (peak5 != null)
            //                                {
            //                                    peak4NeutralMass = Mass.MassFromMz(peak5.X, charge);
            //                                    spacing3 = peak4NeutralMass - peak3NeutralMass;
            //                                    double extendedTheoSpacing = spacingMassRange[2].Mean + spacingMassRange[3].Mean;
            //                                    MassRange doubleSpacing = new MassRange(extendedTheoSpacing, new MassTolerance(MassToleranceType.DA, extendedTheoSpacing * SPACINGPERCENTAGEERROR));

            //                                    // Spacing 1 is bad
            //                                    if (!spacingMassRange[0].Contains(spacing1))
            //                                    {
            //                                        pair.peaks[j, k] = null;
            //                                        peakCount--;

            //                                        // Spacing 1 & 2 are bad
            //                                        if (!spacingMassRange[1].Contains(spacing2))
            //                                        {
            //                                            pair.peaks[j + 1, k] = null;
            //                                            peakCount--;

            //                                            // Spacings 1,2 & 3 are bad
            //                                            if (!doubleSpacing.Contains(spacing3))
            //                                            {
            //                                                pair.peaks[j + 2, k] = null;
            //                                                peakCount--;
            //                                                pair.peaks[j + 4, k] = null;
            //                                                peakCount--;
            //                                            }
            //                                        }
            //                                        // Spacings 1 & 3 are bad
            //                                        else if (!doubleSpacing.Contains(spacing3))
            //                                        {
            //                                            pair.peaks[j + 4, k] = null;
            //                                            peakCount--;
            //                                        }
            //                                    }
            //                                    // Spacing 2 is bad
            //                                    else if (!spacingMassRange[1].Contains(spacing2))
            //                                    {
            //                                        // Spacings 2 & 3 are bad
            //                                        if (!doubleSpacing.Contains(spacing3))
            //                                        {
            //                                            pair.peaks[j + 2, k] = null;
            //                                            peakCount--;
            //                                            pair.peaks[j + 4, k] = null;
            //                                            peakCount--;
            //                                        }
            //                                    }
            //                                    // Spacing 3 is bad
            //                                    else if (!doubleSpacing.Contains(spacing3))
            //                                    {
            //                                        pair.peaks[j + 4, k] = null;
            //                                        peakCount--;
            //                                    }
            //                                }
            //                                // Peaks 1,2,3 & 6
            //                                else if (peak6 != null)
            //                                {
            //                                    peak4NeutralMass = Mass.MassFromMz(peak6.X, charge);
            //                                    spacing3 = peak4NeutralMass - peak3NeutralMass;
            //                                    double extendedTheoSpacing = spacingMassRange[2].Mean + spacingMassRange[3].Mean + spacingMassRange[4].Mean;
            //                                    MassRange tripleSpacing = new MassRange(extendedTheoSpacing, new MassTolerance(MassToleranceType.DA, extendedTheoSpacing * SPACINGPERCENTAGEERROR));

            //                                    // Spacing 1 is bad
            //                                    if (!spacingMassRange[0].Contains(spacing1))
            //                                    {
            //                                        pair.peaks[j, k] = null;
            //                                        peakCount--;

            //                                        // Spacing 1 & 2 are bad
            //                                        if (!spacingMassRange[1].Contains(spacing2))
            //                                        {
            //                                            pair.peaks[j + 1, k] = null;
            //                                            peakCount--;

            //                                            // Spacings 1,2 & 3 are bad
            //                                            if (!tripleSpacing.Contains(spacing3))
            //                                            {
            //                                                pair.peaks[j + 2, k] = null;
            //                                                peakCount--;
            //                                                pair.peaks[j + 5, k] = null;
            //                                                peakCount--;
            //                                            }
            //                                        }
            //                                        // Spacings 1 & 3 are bad
            //                                        else if (!tripleSpacing.Contains(spacing3))
            //                                        {
            //                                            pair.peaks[j + 4, k] = null;
            //                                            peakCount--;
            //                                        }
            //                                    }
            //                                    // Spacing 2 is bad
            //                                    else if (!spacingMassRange[1].Contains(spacing2))
            //                                    {
            //                                        // Spacings 2 & 3 are bad
            //                                        if (!tripleSpacing.Contains(spacing3))
            //                                        {
            //                                            pair.peaks[j + 2, k] = null;
            //                                            peakCount--;
            //                                            pair.peaks[j + 5, k] = null;
            //                                            peakCount--;
            //                                        }
            //                                    }
            //                                    // Spacing 3 is bad
            //                                    else if (!tripleSpacing.Contains(spacing3))
            //                                    {
            //                                        pair.peaks[j + 5, k] = null;
            //                                        peakCount--;
            //                                    }
            //                                }
            //                            }
            //                            else if (peak4 != null)
            //                            {
            //                                peak3NeutralMass = Mass.MassFromMz(peak4.X, charge);
            //                                spacing2 = peak3NeutralMass - peak2NeutralMass;
            //                                double extendedTheoSpacing = spacingMassRange[1].Mean + spacingMassRange[2].Mean;
            //                                MassRange doubleSpacing = new MassRange(extendedTheoSpacing, new MassTolerance(MassToleranceType.DA, extendedTheoSpacing * SPACINGPERCENTAGEERROR));

            //                                // Peaks 1,2,4 & 5
            //                                if (peak5 != null)
            //                                {
            //                                    peak4NeutralMass = Mass.MassFromMz(peak5.X, charge);
            //                                    spacing3 = peak4NeutralMass - peak3NeutralMass;

            //                                    // Spacing 1 is bad
            //                                    if (!spacingMassRange[0].Contains(spacing1))
            //                                    {
            //                                        pair.peaks[j, k] = null;
            //                                        peakCount--;

            //                                        // Spacings 1 & 2 are bad
            //                                        if (!doubleSpacing.Contains(spacing2))
            //                                        {
            //                                            pair.peaks[j + 1, k] = null;
            //                                            peakCount--;

            //                                            // Spacings 1,2 & 3 are bad
            //                                            if (!spacingMassRange[3].Contains(spacing3))
            //                                            {
            //                                                pair.peaks[j + 3, k] = null;
            //                                                peakCount--;
            //                                                pair.peaks[j + 4, k] = null;
            //                                                peakCount--;
            //                                            }
            //                                        }
            //                                        // Spacings 1 & 3 are bad
            //                                        else if (!spacingMassRange[3].Contains(spacing3))
            //                                        {
            //                                            pair.peaks[j + 4, k] = null;
            //                                            peakCount--;
            //                                        }
            //                                    }
            //                                    // Spacing 2 is bad
            //                                    else if (!doubleSpacing.Contains(spacing2))
            //                                    {
            //                                        // Spacings 2 & 3 are bad
            //                                        if (!spacingMassRange[3].Contains(spacing3))
            //                                        {
            //                                            pair.peaks[j + 3, k] = null;
            //                                            peakCount--;
            //                                            pair.peaks[j + 4, k] = null;
            //                                            peakCount--;
            //                                        }
            //                                    }
            //                                    // Spacing 3 is bad
            //                                    else if (!spacingMassRange[3].Contains(spacing3))
            //                                    {
            //                                        pair.peaks[j + 4, k] = null;
            //                                        peakCount--;
            //                                    }
            //                                }
            //                                // Peaks 1,2,4 & 6
            //                                else if (peak6 != null)
            //                                {
            //                                    peak4NeutralMass = Mass.MassFromMz(peak6.X, charge);
            //                                    spacing3 = peak4NeutralMass - peak3NeutralMass;
            //                                    double extendedTheoSpacing1 = spacingMassRange[3].Mean + spacingMassRange[4].Mean;
            //                                    MassRange doubleSpacing1 = new MassRange(extendedTheoSpacing1, new MassTolerance(MassToleranceType.DA, extendedTheoSpacing1 * SPACINGPERCENTAGEERROR));

            //                                    // Spacing 1 is bad
            //                                    if (!spacingMassRange[0].Contains(spacing1))
            //                                    {
            //                                        pair.peaks[j, k] = null;
            //                                        peakCount--;

            //                                        // Spacings 1 & 2 are bad
            //                                        if (!doubleSpacing.Contains(spacing2))
            //                                        {
            //                                            pair.peaks[j + 1, k] = null;
            //                                            peakCount--;

            //                                            // Spacings 1,2 & 3 are bad
            //                                            if (!doubleSpacing1.Contains(spacing3))
            //                                            {
            //                                                pair.peaks[j + 3, k] = null;
            //                                                peakCount--;
            //                                                pair.peaks[j + 5, k] = null;
            //                                                peakCount--;
            //                                            }
            //                                        }
            //                                        // Spacings 1 & 3 are bad
            //                                        else if (!doubleSpacing1.Contains(spacing3))
            //                                        {
            //                                            pair.peaks[j + 5, k] = null;
            //                                            peakCount--;
            //                                        }
            //                                    }
            //                                    // Spacing 2 is bad
            //                                    else if (!doubleSpacing.Contains(spacing2))
            //                                    {
            //                                        // Spacings 2 & 3 are bad
            //                                        if (!doubleSpacing1.Contains(spacing3))
            //                                        {
            //                                            pair.peaks[j + 3, k] = null;
            //                                            peakCount--;
            //                                            pair.peaks[j + 5, k] = null;
            //                                            peakCount--;
            //                                        }
            //                                    }
            //                                    // Spacing 3 is bad
            //                                    else if (!doubleSpacing1.Contains(spacing3))
            //                                    {
            //                                        pair.peaks[j + 5, k] = null;
            //                                        peakCount--;
            //                                    }
            //                                }
            //                            }
            //                            // Peaks 1,2,5 & 6
            //                            else if (peak5 != null)
            //                            {
            //                                peak3NeutralMass = Mass.MassFromMz(peak5.X, charge);
            //                                peak4NeutralMass = Mass.MassFromMz(peak6.X, charge);
            //                                spacing2 = peak3NeutralMass - peak2NeutralMass;
            //                                spacing3 = peak4NeutralMass - peak3NeutralMass;
            //                                double extendedTheoSpacing = spacingMassRange[1].Mean + spacingMassRange[2].Mean + spacingMassRange[3].Mean;
            //                                MassRange tripleSpacing = new MassRange(extendedTheoSpacing, new MassTolerance(MassToleranceType.DA, extendedTheoSpacing * SPACINGPERCENTAGEERROR));

            //                                // Spacing 1 is bad
            //                                if (!spacingMassRange[0].Contains(spacing1))
            //                                {
            //                                    pair.peaks[j, k] = null;
            //                                    peakCount--;

            //                                    // Spacings 1 & 2 are bad
            //                                    if (!tripleSpacing.Contains(spacing2))
            //                                    {
            //                                        pair.peaks[j + 1, k] = null;
            //                                        peakCount--;

            //                                        // Spacings 1,2 & 3 are bad
            //                                        if (!spacingMassRange[4].Contains(spacing3))
            //                                        {
            //                                            pair.peaks[j + 4, k] = null;
            //                                            peakCount--;
            //                                            pair.peaks[j + 5, k] = null;
            //                                            peakCount--;
            //                                        }
            //                                    }
            //                                    // Spacings 1 & 3 are bad
            //                                    else if (!spacingMassRange[4].Contains(spacing3))
            //                                    {
            //                                        pair.peaks[j + 5, k] = null;
            //                                        peakCount--;
            //                                    }
            //                                }
            //                                // Spacing 2 is bad
            //                                else if (!tripleSpacing.Contains(spacing2))
            //                                {
            //                                    // Spacings 2 & 3 are bad
            //                                    if (!spacingMassRange[4].Contains(spacing3))
            //                                    {
            //                                        pair.peaks[j + 4, k] = null;
            //                                        peakCount--;
            //                                        pair.peaks[j + 5, k] = null;
            //                                        peakCount--;
            //                                    }
            //                                }
            //                                // Spacing 3 is bad
            //                                else if (!spacingMassRange[4].Contains(spacing3))
            //                                {
            //                                    pair.peaks[j + 5, k] = null;
            //                                    peakCount--;
            //                                }
            //                            }
            //                        }
            //                        else if (peak3 != null)
            //                        {
            //                            peak2NeutralMass = Mass.MassFromMz(peak3.X, charge);
            //                            spacing1 = peak2NeutralMass - peak1NeutralMass;
            //                            double extendedTheoSpacing = spacingMassRange[0].Mean + spacingMassRange[1].Mean;
            //                            MassRange doubleSpacing = new MassRange(extendedTheoSpacing, new MassTolerance(MassToleranceType.DA, extendedTheoSpacing * SPACINGPERCENTAGEERROR));

            //                            if (peak4 != null)
            //                            {
            //                                peak3NeutralMass = Mass.MassFromMz(peak4.X, charge);
            //                                spacing2 = peak3NeutralMass - peak2NeutralMass;

            //                                // Peaks 1,3,4 & 5
            //                                if (peak5 != null)
            //                                {
            //                                    peak4NeutralMass = Mass.MassFromMz(peak5.X, charge);
            //                                    spacing3 = peak4NeutralMass - peak3NeutralMass;

            //                                    // Spacing 1 is bad
            //                                    if (!doubleSpacing.Contains(spacing1))
            //                                    {
            //                                        pair.peaks[j, k] = null;
            //                                        peakCount--;

            //                                        // Spacings 1 & 2 are bad
            //                                        if (!spacingMassRange[2].Contains(spacing2))
            //                                        {
            //                                            pair.peaks[j + 2, k] = null;
            //                                            peakCount--;

            //                                            // Spacings 1,2 & 3 are bad
            //                                            if (!spacingMassRange[3].Contains(spacing3))
            //                                            {
            //                                                pair.peaks[j + 3, k] = null;
            //                                                peakCount--;
            //                                                pair.peaks[j + 4, k] = null;
            //                                                peakCount--;
            //                                            }
            //                                        }
            //                                        // Spacings 1 & 3 are bad
            //                                        else if (!spacingMassRange[3].Contains(spacing3))
            //                                        {
            //                                            pair.peaks[j + 4, k] = null;
            //                                            peakCount--;
            //                                        }
            //                                    }
            //                                    // Spacing 2 is bad
            //                                    else if (!spacingMassRange[2].Contains(spacing2))
            //                                    {
            //                                        // Spacings 2 & 3 are bad
            //                                        if (!spacingMassRange[3].Contains(spacing3))
            //                                        {
            //                                            pair.peaks[j + 3, k] = null;
            //                                            peakCount--;
            //                                            pair.peaks[j + 4, k] = null;
            //                                            peakCount--;
            //                                        }
            //                                    }
            //                                    // Spacing 3 is bad
            //                                    else if (!spacingMassRange[3].Contains(spacing3))
            //                                    {
            //                                        pair.peaks[j + 4, k] = null;
            //                                        peakCount--;
            //                                    }
            //                                }
            //                                // Peaks 1,3,4 & 6
            //                                else if (peak6 != null)
            //                                {
            //                                    peak4NeutralMass = Mass.MassFromMz(peak6.X, charge);
            //                                    spacing3 = peak4NeutralMass - peak3NeutralMass;
            //                                    double extendedTheoSpacing1 = spacingMassRange[3].Mean + spacingMassRange[4].Mean;
            //                                    MassRange doubleSpacing1 = new MassRange(extendedTheoSpacing1, new MassTolerance(MassToleranceType.DA, extendedTheoSpacing1 * SPACINGPERCENTAGEERROR));

            //                                    // Spacing 1 is bad
            //                                    if (!doubleSpacing.Contains(spacing1))
            //                                    {
            //                                        pair.peaks[j, k] = null;
            //                                        peakCount--;

            //                                        // Spacings 1 & 2 are bad
            //                                        if (!spacingMassRange[2].Contains(spacing2))
            //                                        {
            //                                            pair.peaks[j + 2, k] = null;
            //                                            peakCount--;

            //                                            // Spacings 1,2 & 3 are bad
            //                                            if (!doubleSpacing1.Contains(spacing3))
            //                                            {
            //                                                pair.peaks[j + 3, k] = null;
            //                                                peakCount--;
            //                                                pair.peaks[j + 5, k] = null;
            //                                                peakCount--;
            //                                            }
            //                                        }
            //                                        // Spacings 1 & 3 are bad
            //                                        else if (!doubleSpacing1.Contains(spacing3))
            //                                        {
            //                                            pair.peaks[j + 5, k] = null;
            //                                            peakCount--;
            //                                        }
            //                                    }
            //                                    // Spacing 2 is bad
            //                                    else if (!spacingMassRange[2].Contains(spacing2))
            //                                    {
            //                                        // Spacings 2 & 3 are bad
            //                                        if (!doubleSpacing1.Contains(spacing3))
            //                                        {
            //                                            pair.peaks[j + 3, k] = null;
            //                                            peakCount--;
            //                                            pair.peaks[j + 5, k] = null;
            //                                            peakCount--;
            //                                        }
            //                                    }
            //                                    // Spacing 3 is bad
            //                                    else if (!doubleSpacing1.Contains(spacing3))
            //                                    {
            //                                        pair.peaks[j + 5, k] = null;
            //                                        peakCount--;
            //                                    }
            //                                }
            //                            }
            //                            // Peaks 1,3,5 & 6
            //                            else if (peak5 != null)
            //                            {
            //                                peak3NeutralMass = Mass.MassFromMz(peak5.X, charge);
            //                                peak4NeutralMass = Mass.MassFromMz(peak6.X, charge);
            //                                spacing2 = peak3NeutralMass - peak2NeutralMass;
            //                                spacing3 = peak4NeutralMass - peak3NeutralMass;
            //                                double extendedTheoSpacing1 = spacingMassRange[2].Mean + spacingMassRange[3].Mean;
            //                                MassRange doubleSpacing1 = new MassRange(extendedTheoSpacing1, new MassTolerance(MassToleranceType.DA, extendedTheoSpacing1 * SPACINGPERCENTAGEERROR));

            //                                // Spacing 1 is bad
            //                                if (!doubleSpacing.Contains(spacing1))
            //                                {
            //                                    pair.peaks[j, k] = null;
            //                                    peakCount--;

            //                                    // Spacings 1 & 2 are bad
            //                                    if (!doubleSpacing1.Contains(spacing2))
            //                                    {
            //                                        pair.peaks[j + 2, k] = null;
            //                                        peakCount--;

            //                                        // Spacings 1,2 & 3 are bad
            //                                        if (!spacingMassRange[4].Contains(spacing3))
            //                                        {
            //                                            pair.peaks[j + 4, k] = null;
            //                                            peakCount--;
            //                                            pair.peaks[j + 5, k] = null;
            //                                            peakCount--;
            //                                        }
            //                                    }
            //                                    // Spacings 1 & 3 are bad
            //                                    else if (!spacingMassRange[4].Contains(spacing3))
            //                                    {
            //                                        pair.peaks[j + 5, k] = null;
            //                                        peakCount--;
            //                                    }
            //                                }
            //                                // Spacing 2 is bad
            //                                else if (!doubleSpacing1.Contains(spacing2))
            //                                {
            //                                    // Spacings 2 & 3 are bad
            //                                    if (!spacingMassRange[4].Contains(spacing3))
            //                                    {
            //                                        pair.peaks[j + 4, k] = null;
            //                                        peakCount--;
            //                                        pair.peaks[j + 5, k] = null;
            //                                        peakCount--;
            //                                    }
            //                                }
            //                                // Spacing 3 is bad
            //                                else if (!spacingMassRange[4].Contains(spacing3))
            //                                {
            //                                    pair.peaks[j + 5, k] = null;
            //                                    peakCount--;
            //                                }
            //                            }
            //                        }
            //                        // Peaks 1,4,5 & 6
            //                        else if (peak4 != null)
            //                        {
            //                            peak2NeutralMass = Mass.MassFromMz(peak4.X, charge);
            //                            peak3NeutralMass = Mass.MassFromMz(peak5.X, charge);
            //                            peak4NeutralMass = Mass.MassFromMz(peak6.X, charge);
            //                            spacing1 = peak2NeutralMass - peak1NeutralMass;
            //                            spacing2 = peak3NeutralMass - peak2NeutralMass;
            //                            spacing3 = peak4NeutralMass - peak3NeutralMass;
            //                            double extendedTheoSpacing = spacingMassRange[1].Mean + spacingMassRange[2].Mean + spacingMassRange[3].Mean;
            //                            MassRange tripleSpacing = new MassRange(extendedTheoSpacing, new MassTolerance(MassToleranceType.DA, extendedTheoSpacing * SPACINGPERCENTAGEERROR));

            //                            // Spacing 1 is bad
            //                            if (!tripleSpacing.Contains(spacing1))
            //                            {
            //                                pair.peaks[j, k] = null;
            //                                peakCount--;

            //                                // Spacings 1 & 2 are bad
            //                                if (!spacingMassRange[3].Contains(spacing2))
            //                                {
            //                                    pair.peaks[j + 3, k] = null;
            //                                    peakCount--;

            //                                    // Spacings 1,2 & 3 are bad
            //                                    if (!spacingMassRange[4].Contains(spacing3))
            //                                    {
            //                                        pair.peaks[j + 4, k] = null;
            //                                        peakCount--;
            //                                        pair.peaks[j + 5, k] = null;
            //                                        peakCount--;
            //                                    }
            //                                }
            //                                // Spacings 1 & 3 are bad
            //                                else if (!spacingMassRange[4].Contains(spacing3))
            //                                {
            //                                    pair.peaks[j + 5, k] = null;
            //                                    peakCount--;
            //                                }
            //                            }
            //                            // Spacing 2 is bad
            //                            else if (!spacingMassRange[3].Contains(spacing2))
            //                            {
            //                                // Spacings 2 & 3 are bad
            //                                if (!spacingMassRange[4].Contains(spacing3))
            //                                {
            //                                    pair.peaks[j + 4, k] = null;
            //                                    peakCount--;
            //                                    pair.peaks[j + 5, k] = null;
            //                                    peakCount--;
            //                                }
            //                            }
            //                            // Spacing 3 is bad
            //                            else if (!spacingMassRange[4].Contains(spacing3))
            //                            {
            //                                pair.peaks[j + 5, k] = null;
            //                                peakCount--;
            //                            }
            //                        }
            //                    }
            //                    else if (peak2 != null)
            //                    {
            //                        peak1NeutralMass = Mass.MassFromMz(peak2.X, charge);

            //                        if (peak3 != null)
            //                        {
            //                            peak2NeutralMass = Mass.MassFromMz(peak3.X, charge);
            //                            spacing1 = peak2NeutralMass - peak1NeutralMass;

            //                            if (peak4 != null)
            //                            {
            //                                peak3NeutralMass = Mass.MassFromMz(peak4.X, charge);
            //                                spacing2 = peak3NeutralMass - peak2NeutralMass;

            //                                // Peaks 2,3,4 & 5
            //                                if (peak5 != null)
            //                                {
            //                                    peak4NeutralMass = Mass.MassFromMz(peak5.X, charge);
            //                                    spacing3 = peak4NeutralMass - peak3NeutralMass;

            //                                    // Spacing 1 is bad
            //                                    if (!spacingMassRange[1].Contains(spacing1))
            //                                    {
            //                                        pair.peaks[j + 1, k] = null;
            //                                        peakCount--;

            //                                        // Spacings 1 & 2 are bad
            //                                        if (!spacingMassRange[2].Contains(spacing2))
            //                                        {
            //                                            pair.peaks[j + 2, k] = null;
            //                                            peakCount--;

            //                                            // Spacings 1,2 & 3 are bad
            //                                            if (!spacingMassRange[3].Contains(spacing3))
            //                                            {
            //                                                pair.peaks[j + 3, k] = null;
            //                                                peakCount--;
            //                                                pair.peaks[j + 4, k] = null;
            //                                                peakCount--;
            //                                            }
            //                                        }
            //                                        // Spacings 1 & 3 are bad
            //                                        else if (!spacingMassRange[3].Contains(spacing3))
            //                                        {
            //                                            pair.peaks[j + 4, k] = null;
            //                                            peakCount--;
            //                                        }
            //                                    }
            //                                    // Spacing 2 is bad
            //                                    else if (!spacingMassRange[2].Contains(spacing2))
            //                                    {
            //                                        // Spacings 2 & 3 are bad
            //                                        if (!spacingMassRange[3].Contains(spacing3))
            //                                        {
            //                                            pair.peaks[j + 3, k] = null;
            //                                            peakCount--;
            //                                            pair.peaks[j + 4, k] = null;
            //                                            peakCount--;
            //                                        }
            //                                    }
            //                                    // Spacing 3 is bad
            //                                    else if (!spacingMassRange[3].Contains(spacing3))
            //                                    {
            //                                        pair.peaks[j + 4, k] = null;
            //                                        peakCount--;
            //                                    }
            //                                }
            //                                // Peaks 2,3,4 & 6
            //                                else if (peak6 != null)
            //                                {
            //                                    peak4NeutralMass = Mass.MassFromMz(peak6.X, charge);
            //                                    spacing3 = peak4NeutralMass - peak3NeutralMass;
            //                                    double extendedTheoSpacing = spacingMassRange[3].Mean + spacingMassRange[4].Mean;
            //                                    MassRange doubleSpacing = new MassRange(extendedTheoSpacing, new MassTolerance(MassToleranceType.DA, extendedTheoSpacing * SPACINGPERCENTAGEERROR));

            //                                    // Spacing 1 is bad
            //                                    if (!spacingMassRange[1].Contains(spacing1))
            //                                    {
            //                                        pair.peaks[j + 1, k] = null;
            //                                        peakCount--;

            //                                        // Spacings 1 & 2 are bad
            //                                        if (!spacingMassRange[2].Contains(spacing2))
            //                                        {
            //                                            pair.peaks[j + 2, k] = null;
            //                                            peakCount--;

            //                                            // Spacings 1,2 & 3 are bad
            //                                            if (!doubleSpacing.Contains(spacing3))
            //                                            {
            //                                                pair.peaks[j + 3, k] = null;
            //                                                peakCount--;
            //                                                pair.peaks[j + 5, k] = null;
            //                                                peakCount--;
            //                                            }
            //                                        }
            //                                        // Spacings 1 & 3 are bad
            //                                        else if (!doubleSpacing.Contains(spacing3))
            //                                        {
            //                                            pair.peaks[j + 5, k] = null;
            //                                            peakCount--;
            //                                        }
            //                                    }
            //                                    // Spacing 2 is bad
            //                                    else if (!spacingMassRange[2].Contains(spacing2))
            //                                    {
            //                                        // Spacings 2 & 3 are bad
            //                                        if (!doubleSpacing.Contains(spacing3))
            //                                        {
            //                                            pair.peaks[j + 3, k] = null;
            //                                            peakCount--;
            //                                            pair.peaks[j + 5, k] = null;
            //                                            peakCount--;
            //                                        }
            //                                    }
            //                                    // Spacing 3 is bad
            //                                    else if (!doubleSpacing.Contains(spacing3))
            //                                    {
            //                                        pair.peaks[j + 5, k] = null;
            //                                        peakCount--;
            //                                    }
            //                                }
            //                            }
            //                            // Peaks 2,3,5 & 6
            //                            else if (peak5 != null)
            //                            {
            //                                peak3NeutralMass = Mass.MassFromMz(peak5.X, charge);
            //                                peak4NeutralMass = Mass.MassFromMz(peak6.X, charge);
            //                                spacing2 = peak3NeutralMass - peak2NeutralMass;
            //                                spacing3 = peak4NeutralMass - peak3NeutralMass;
            //                                double extendedTheoSpacing = spacingMassRange[2].Mean + spacingMassRange[3].Mean;
            //                                MassRange doubleSpacing = new MassRange(extendedTheoSpacing, new MassTolerance(MassToleranceType.DA, extendedTheoSpacing * SPACINGPERCENTAGEERROR));

            //                                // Spacing 1 is bad
            //                                if (!spacingMassRange[1].Contains(spacing1))
            //                                {
            //                                    pair.peaks[j + 1, k] = null;
            //                                    peakCount--;

            //                                    // Spacings 1 & 2 are bad
            //                                    if (!doubleSpacing.Contains(spacing2))
            //                                    {
            //                                        pair.peaks[j + 2, k] = null;
            //                                        peakCount--;

            //                                        // Spacings 1,2 & 3 are bad
            //                                        if (!spacingMassRange[4].Contains(spacing3))
            //                                        {
            //                                            pair.peaks[j + 4, k] = null;
            //                                            peakCount--;
            //                                            pair.peaks[j + 5, k] = null;
            //                                            peakCount--;
            //                                        }
            //                                    }
            //                                    // Spacings 1 & 3 are bad
            //                                    else if (!spacingMassRange[4].Contains(spacing3))
            //                                    {
            //                                        pair.peaks[j + 5, k] = null;
            //                                        peakCount--;
            //                                    }
            //                                }
            //                                // Spacing 2 is bad
            //                                else if (!doubleSpacing.Contains(spacing2))
            //                                {
            //                                    // Spacings 2 & 3 are bad
            //                                    if (!spacingMassRange[4].Contains(spacing3))
            //                                    {
            //                                        pair.peaks[j + 4, k] = null;
            //                                        peakCount--;
            //                                        pair.peaks[j + 5, k] = null;
            //                                        peakCount--;
            //                                    }
            //                                }
            //                                // Spacing 3 is bad
            //                                else if (!spacingMassRange[4].Contains(spacing3))
            //                                {
            //                                    pair.peaks[j + 5, k] = null;
            //                                    peakCount--;
            //                                }
            //                            }
            //                        }
            //                        // Peaks 2,4,5 & 6
            //                        else if (peak4 != null)
            //                        {
            //                            peak2NeutralMass = Mass.MassFromMz(peak4.X, charge);
            //                            peak3NeutralMass = Mass.MassFromMz(peak5.X, charge);
            //                            peak4NeutralMass = Mass.MassFromMz(peak6.X, charge);
            //                            spacing1 = peak2NeutralMass - peak1NeutralMass;
            //                            spacing2 = peak3NeutralMass - peak2NeutralMass;
            //                            spacing3 = peak4NeutralMass - peak3NeutralMass;
            //                            double extendedTheoSpacing = spacingMassRange[1].Mean + spacingMassRange[2].Mean;
            //                            MassRange doubleSpacing = new MassRange(extendedTheoSpacing, new MassTolerance(MassToleranceType.DA, extendedTheoSpacing * SPACINGPERCENTAGEERROR));

            //                            // Spacing 1 is bad
            //                            if (!doubleSpacing.Contains(spacing1))
            //                            {
            //                                pair.peaks[j + 1, k] = null;
            //                                peakCount--;

            //                                // Spacings 1 & 2 are bad
            //                                if (!spacingMassRange[3].Contains(spacing2))
            //                                {
            //                                    pair.peaks[j + 3, k] = null;
            //                                    peakCount--;

            //                                    // Spacings 1,2 & 3 are bad
            //                                    if (!spacingMassRange[4].Contains(spacing3))
            //                                    {
            //                                        pair.peaks[j + 4, k] = null;
            //                                        peakCount--;
            //                                        pair.peaks[j + 5, k] = null;
            //                                        peakCount--;
            //                                    }
            //                                }
            //                                // Spacings 1 & 3 are bad
            //                                else if (!spacingMassRange[4].Contains(spacing3))
            //                                {
            //                                    pair.peaks[j + 5, k] = null;
            //                                    peakCount--;
            //                                }
            //                            }
            //                            // Spacing 2 is bad
            //                            else if (!spacingMassRange[3].Contains(spacing2))
            //                            {
            //                                // Spacings 2 & 3 are bad
            //                                if (!spacingMassRange[4].Contains(spacing3))
            //                                {
            //                                    pair.peaks[j + 4, k] = null;
            //                                    peakCount--;
            //                                    pair.peaks[j + 5, k] = null;
            //                                    peakCount--;
            //                                }
            //                            }
            //                            // Spacing 3 is bad
            //                            else if (!spacingMassRange[4].Contains(spacing3))
            //                            {
            //                                pair.peaks[j + 5, k] = null;
            //                                peakCount--;
            //                            }
            //                        }
            //                    }
            //                    // Peaks 3,4,5 & 6
            //                    else if (peak3 != null)
            //                    {
            //                        peak1NeutralMass = Mass.MassFromMz(peak3.X, charge);
            //                        peak2NeutralMass = Mass.MassFromMz(peak4.X, charge);
            //                        peak3NeutralMass = Mass.MassFromMz(peak5.X, charge);
            //                        peak4NeutralMass = Mass.MassFromMz(peak6.X, charge);
            //                        spacing1 = peak2NeutralMass - peak1NeutralMass;
            //                        spacing2 = peak3NeutralMass - peak2NeutralMass;
            //                        spacing3 = peak4NeutralMass - peak3NeutralMass;

            //                        // Spacing 1 is bad
            //                        if (!spacingMassRange[2].Contains(spacing1))
            //                        {
            //                            pair.peaks[j + 2, k] = null;
            //                            peakCount--;

            //                            // Spacings 1 & 2 are bad
            //                            if (!spacingMassRange[3].Contains(spacing2))
            //                            {
            //                                pair.peaks[j + 3, k] = null;
            //                                peakCount--;

            //                                // Spacings 1,2 & 3 are bad
            //                                if (!spacingMassRange[4].Contains(spacing3))
            //                                {
            //                                    pair.peaks[j + 4, k] = null;
            //                                    peakCount--;
            //                                    pair.peaks[j + 5, k] = null;
            //                                    peakCount--;
            //                                }
            //                            }
            //                            // Spacings 1 & 3 are bad
            //                            else if (!spacingMassRange[4].Contains(spacing3))
            //                            {
            //                                pair.peaks[j + 5, k] = null;
            //                                peakCount--;
            //                            }
            //                        }
            //                        // Spacing 2 is bad
            //                        else if (!spacingMassRange[3].Contains(spacing2))
            //                        {
            //                            // Spacings 2 & 3 are bad
            //                            if (!spacingMassRange[4].Contains(spacing3))
            //                            {
            //                                pair.peaks[j + 4, k] = null;
            //                                peakCount--;
            //                                pair.peaks[j + 5, k] = null;
            //                                peakCount--;
            //                            }
            //                        }
            //                        // Spacing 3 is bad
            //                        else if (!spacingMassRange[4].Contains(spacing3))
            //                        {
            //                            pair.peaks[j + 5, k] = null;
            //                            peakCount--;
            //                        }
            //                    }
            //                }
            //                else if (Form1.NOISEBANDCAP && peakCount == 3)
            //                {
            //                    if (peak1 != null)
            //                    {
            //                        peak1NeutralMass = Mass.MassFromMz(peak1.X, charge);

            //                        if (peak2 != null)
            //                        {
            //                            peak2NeutralMass = Mass.MassFromMz(peak2.X, charge);
            //                            spacing1 = peak2NeutralMass - peak1NeutralMass;

            //                            // Peaks 1,2 & 3
            //                            if (peak3 != null)
            //                            {
            //                                peak3NeutralMass = Mass.MassFromMz(peak3.X, charge);
            //                                spacing2 = peak3NeutralMass - peak2NeutralMass;

            //                                // Spacing 1 is bad
            //                                if (!spacingMassRange[0].Contains(spacing1))
            //                                {
            //                                    pair.peaks[j, k] = null;
            //                                    peakCount--;

            //                                    // Spacing 1 & 2 are bad
            //                                    if (!spacingMassRange[1].Contains(spacing2))
            //                                    {
            //                                        pair.peaks[j + 1, k] = null;
            //                                        peakCount--;
            //                                        pair.peaks[j + 2, k] = null;
            //                                        peakCount--;
            //                                    }
            //                                }
            //                                // Spacing 2 is bad
            //                                else if (!spacingMassRange[1].Contains(spacing2))
            //                                {
            //                                    pair.peaks[j + 2, k] = null;
            //                                    peakCount--;
            //                                }
            //                            }
            //                            // Peaks 1,2 & 4
            //                            else if (peak4 != null)
            //                            {
            //                                peak3NeutralMass = Mass.MassFromMz(peak4.X, charge);
            //                                spacing2 = peak3NeutralMass - peak2NeutralMass;
            //                                double extendedTheoSpacing = spacingMassRange[1].Mean + spacingMassRange[2].Mean;
            //                                MassRange doubleSpacing = new MassRange(extendedTheoSpacing, new MassTolerance(MassToleranceType.DA, extendedTheoSpacing * SPACINGPERCENTAGEERROR));

            //                                // Spacing 1 is bad
            //                                if (!spacingMassRange[0].Contains(spacing1))
            //                                {
            //                                    pair.peaks[j, k] = null;
            //                                    peakCount--;

            //                                    // Spacing 1 & 2 are bad
            //                                    if (!doubleSpacing.Contains(spacing2))
            //                                    {
            //                                        pair.peaks[j + 1, k] = null;
            //                                        peakCount--;
            //                                        pair.peaks[j + 3, k] = null;
            //                                        peakCount--;
            //                                    }
            //                                }
            //                                // Spacing 2 is bad
            //                                else if (!doubleSpacing.Contains(spacing2))
            //                                {
            //                                    pair.peaks[j + 3, k] = null;
            //                                    peakCount--;
            //                                }
            //                            }
            //                            // Peaks 1,2 & 5
            //                            else if (peak5 != null)
            //                            {
            //                                peak3NeutralMass = Mass.MassFromMz(peak5.X, charge);
            //                                spacing2 = peak3NeutralMass - peak2NeutralMass;
            //                                double extendedTheoSpacing = spacingMassRange[1].Mean + spacingMassRange[2].Mean + spacingMassRange[3].Mean;
            //                                MassRange tripleSpacing = new MassRange(extendedTheoSpacing, new MassTolerance(MassToleranceType.DA, extendedTheoSpacing * SPACINGPERCENTAGEERROR));

            //                                // Spacing 1 is bad
            //                                if (!spacingMassRange[0].Contains(spacing1))
            //                                {
            //                                    pair.peaks[j, k] = null;
            //                                    peakCount--;

            //                                    // Spacing 1 & 2 are bad
            //                                    if (!tripleSpacing.Contains(spacing2))
            //                                    {
            //                                        pair.peaks[j + 1, k] = null;
            //                                        peakCount--;
            //                                        pair.peaks[j + 4, k] = null;
            //                                        peakCount--;
            //                                    }
            //                                }
            //                                // Spacing 2 is bad
            //                                else if (!tripleSpacing.Contains(spacing2))
            //                                {
            //                                    pair.peaks[j + 4, k] = null;
            //                                    peakCount--;
            //                                }
            //                            }
            //                            // Peaks 1,2 & 6
            //                            else if (peak6 != null)
            //                            {
            //                                peak3NeutralMass = Mass.MassFromMz(peak6.X, charge);
            //                                spacing2 = peak3NeutralMass - peak2NeutralMass;
            //                                double extendedTheoSpacing = spacingMassRange[1].Mean + spacingMassRange[2].Mean + spacingMassRange[3].Mean + spacingMassRange[4].Mean;
            //                                MassRange quadrupleSpacing = new MassRange(extendedTheoSpacing, new MassTolerance(MassToleranceType.DA, extendedTheoSpacing * SPACINGPERCENTAGEERROR));

            //                                // Spacing 1 is bad
            //                                if (!spacingMassRange[0].Contains(spacing1))
            //                                {
            //                                    pair.peaks[j, k] = null;
            //                                    peakCount--;

            //                                    // Spacing 1 & 2 are bad
            //                                    if (!quadrupleSpacing.Contains(spacing2))
            //                                    {
            //                                        pair.peaks[j + 1, k] = null;
            //                                        peakCount--;
            //                                        pair.peaks[j + 5, k] = null;
            //                                        peakCount--;
            //                                    }
            //                                }
            //                                // Spacing 2 is bad
            //                                else if (!quadrupleSpacing.Contains(spacing2))
            //                                {
            //                                    pair.peaks[j + 5, k] = null;
            //                                    peakCount--;
            //                                }
            //                            }
            //                        }
            //                        else if (peak3 != null)
            //                        {
            //                            peak2NeutralMass = Mass.MassFromMz(peak3.X, charge);
            //                            spacing1 = peak2NeutralMass - peak1NeutralMass;
            //                            double extendedTheoSpacing = spacingMassRange[0].Mean + spacingMassRange[1].Mean;
            //                            MassRange doubleSpacing = new MassRange(extendedTheoSpacing, new MassTolerance(MassToleranceType.DA, extendedTheoSpacing * SPACINGPERCENTAGEERROR));

            //                            // Peaks 1,3 & 4
            //                            if (peak4 != null)
            //                            {
            //                                peak3NeutralMass = Mass.MassFromMz(peak4.X, charge);
            //                                spacing2 = peak3NeutralMass - peak2NeutralMass;

            //                                // Spacing 1 is bad
            //                                if (!doubleSpacing.Contains(spacing1))
            //                                {
            //                                    pair.peaks[j, k] = null;
            //                                    peakCount--;

            //                                    // Spacing 1 & 2 are bad
            //                                    if (!spacingMassRange[2].Contains(spacing2))
            //                                    {
            //                                        pair.peaks[j + 2, k] = null;
            //                                        peakCount--;
            //                                        pair.peaks[j + 3, k] = null;
            //                                        peakCount--;
            //                                    }
            //                                }
            //                                // Spacing 2 is bad
            //                                else if (!spacingMassRange[2].Contains(spacing2))
            //                                {
            //                                    pair.peaks[j + 3, k] = null;
            //                                    peakCount--;
            //                                }
            //                            }
            //                            // Peaks 1,3 & 5
            //                            else if (peak5 != null)
            //                            {
            //                                peak3NeutralMass = Mass.MassFromMz(peak5.X, charge);
            //                                spacing2 = peak3NeutralMass - peak2NeutralMass;
            //                                double extendedTheoSpacing1 = spacingMassRange[2].Mean + spacingMassRange[3].Mean;
            //                                MassRange doubleSpacing1 = new MassRange(extendedTheoSpacing1, new MassTolerance(MassToleranceType.DA, extendedTheoSpacing1 * SPACINGPERCENTAGEERROR));

            //                                // Spacing 1 is bad
            //                                if (!doubleSpacing.Contains(spacing1))
            //                                {
            //                                    pair.peaks[j, k] = null;
            //                                    peakCount--;

            //                                    // Spacing 1 & 2 are bad
            //                                    if (!doubleSpacing1.Contains(spacing2))
            //                                    {
            //                                        pair.peaks[j + 2, k] = null;
            //                                        peakCount--;
            //                                        pair.peaks[j + 4, k] = null;
            //                                        peakCount--;
            //                                    }
            //                                }
            //                                // Spacing 2 is bad
            //                                else if (!doubleSpacing1.Contains(spacing2))
            //                                {
            //                                    pair.peaks[j + 4, k] = null;
            //                                    peakCount--;
            //                                }
            //                            }
            //                            // Peaks 1,3 & 6
            //                            else if (peak6 != null)
            //                            {
            //                                peak3NeutralMass = Mass.MassFromMz(peak6.X, charge);
            //                                spacing2 = peak3NeutralMass - peak2NeutralMass;
            //                                double extendedTheoSpacing2 = spacingMassRange[2].Mean + spacingMassRange[3].Mean + spacingMassRange[4].Mean;
            //                                MassRange tripleSpacing = new MassRange(extendedTheoSpacing2, new MassTolerance(MassToleranceType.DA, extendedTheoSpacing2 * SPACINGPERCENTAGEERROR));

            //                                // Spacing 1 is bad
            //                                if (!doubleSpacing.Contains(spacing1))
            //                                {
            //                                    pair.peaks[j, k] = null;
            //                                    peakCount--;

            //                                    // Spacing 1 & 2 are bad
            //                                    if (!tripleSpacing.Contains(spacing2))
            //                                    {
            //                                        pair.peaks[j + 2, k] = null;
            //                                        peakCount--;
            //                                        pair.peaks[j + 5, k] = null;
            //                                        peakCount--;
            //                                    }
            //                                }
            //                                // Spacing 2 is bad
            //                                else if (!tripleSpacing.Contains(spacing2))
            //                                {
            //                                    pair.peaks[j + 5, k] = null;
            //                                    peakCount--;
            //                                }
            //                            }
            //                        }
            //                        else if (peak4 != null)
            //                        {
            //                            peak2NeutralMass = Mass.MassFromMz(peak4.X, charge);
            //                            spacing1 = peak2NeutralMass - peak1NeutralMass;
            //                            double extendedTheoSpacing = spacingMassRange[0].Mean + spacingMassRange[1].Mean + spacingMassRange[2].Mean;
            //                            MassRange tripleSpacing = new MassRange(extendedTheoSpacing, new MassTolerance(MassToleranceType.DA, extendedTheoSpacing * SPACINGPERCENTAGEERROR));

            //                            // Peaks 1,4 & 5
            //                            if (peak5 != null)
            //                            {
            //                                peak3NeutralMass = Mass.MassFromMz(peak5.X, charge);
            //                                spacing2 = peak3NeutralMass - peak2NeutralMass;

            //                                // Spacing 1 is bad
            //                                if (!tripleSpacing.Contains(spacing1))
            //                                {
            //                                    pair.peaks[j, k] = null;
            //                                    peakCount--;

            //                                    // Spacings 1 & 2 are bad
            //                                    if (!spacingMassRange[3].Contains(spacing2))
            //                                    {
            //                                        pair.peaks[j + 3, k] = null;
            //                                        peakCount--;
            //                                        pair.peaks[j + 4, k] = null;
            //                                        peakCount--;
            //                                    }
            //                                }
            //                                // Spacing 2 is bad
            //                                else if (!spacingMassRange[3].Contains(spacing2))
            //                                {
            //                                    pair.peaks[j + 4, k] = null;
            //                                    peakCount--;
            //                                }
            //                            }
            //                            // Peaks 1,4 & 6
            //                            else if (peak6 != null)
            //                            {
            //                                peak3NeutralMass = Mass.MassFromMz(peak6.X, charge);
            //                                spacing2 = peak3NeutralMass - peak2NeutralMass;
            //                                double extendedTheoSpacing1 = spacingMassRange[3].Mean + spacingMassRange[4].Mean;
            //                                MassRange doubleSpacing = new MassRange(extendedTheoSpacing1, new MassTolerance(MassToleranceType.DA, extendedTheoSpacing1 * SPACINGPERCENTAGEERROR));

            //                                // Spacing 1 is bad
            //                                if (!tripleSpacing.Contains(spacing1))
            //                                {
            //                                    pair.peaks[j, k] = null;
            //                                    peakCount--;

            //                                    // Spacings 1 & 2 are bad
            //                                    if (!doubleSpacing.Contains(spacing2))
            //                                    {
            //                                        pair.peaks[j + 3, k] = null;
            //                                        peakCount--;
            //                                        pair.peaks[j + 5, k] = null;
            //                                        peakCount--;
            //                                    }
            //                                }
            //                                // Spacing 2 is bad
            //                                else if (!doubleSpacing.Contains(spacing2))
            //                                {
            //                                    pair.peaks[j + 5, k] = null;
            //                                    peakCount--;
            //                                }
            //                            }
            //                        }
            //                        // Peaks 1,5 & 6
            //                        else if (peak5 != null)
            //                        {
            //                            peak2NeutralMass = Mass.MassFromMz(peak5.X, charge);
            //                            peak3NeutralMass = Mass.MassFromMz(peak6.X, charge);
            //                            spacing1 = peak2NeutralMass - peak1NeutralMass;
            //                            spacing2 = peak3NeutralMass - peak2NeutralMass;
            //                            double extendedTheoSpacing = spacingMassRange[0].Mean + spacingMassRange[1].Mean + spacingMassRange[2].Mean + spacingMassRange[3].Mean;
            //                            MassRange quadrupleSpacing = new MassRange(extendedTheoSpacing, new MassTolerance(MassToleranceType.DA, extendedTheoSpacing * SPACINGPERCENTAGEERROR));

            //                            // Spacing 1 is bad
            //                            if (!quadrupleSpacing.Contains(spacing1))
            //                            {
            //                                pair.peaks[j, k] = null;
            //                                peakCount--;

            //                                // Spacings 1 & 2 are bad
            //                                if (!spacingMassRange[4].Contains(spacing2))
            //                                {
            //                                    pair.peaks[j + 4, k] = null;
            //                                    peakCount--;
            //                                    pair.peaks[j + 5, k] = null;
            //                                    peakCount--;
            //                                }
            //                            }
            //                            // Spacing 2 is bad
            //                            else if (!spacingMassRange[4].Contains(spacing2))
            //                            {
            //                                pair.peaks[j + 5, k] = null;
            //                                peakCount--;
            //                            }
            //                        }
            //                    }
            //                    else if (peak2 != null)
            //                    {
            //                        peak1NeutralMass = Mass.MassFromMz(peak2.X, charge);

            //                        if (peak3 != null)
            //                        {
            //                            peak2NeutralMass = Mass.MassFromMz(peak3.X, charge);
            //                            spacing1 = peak2NeutralMass - peak1NeutralMass;

            //                            // Peaks 2,3 & 4
            //                            if (peak4 != null)
            //                            {
            //                                peak3NeutralMass = Mass.MassFromMz(peak4.X, charge);
            //                                spacing2 = peak3NeutralMass - peak2NeutralMass;

            //                                // Spacing 1 is bad
            //                                if (!spacingMassRange[1].Contains(spacing1))
            //                                {
            //                                    pair.peaks[j + 1, k] = null;
            //                                    peakCount--;

            //                                    // Spacings 1 & 2 are bad
            //                                    if (!spacingMassRange[2].Contains(spacing2))
            //                                    {
            //                                        pair.peaks[j + 2, k] = null;
            //                                        peakCount--;
            //                                        pair.peaks[j + 3, k] = null;
            //                                        peakCount--;
            //                                    }
            //                                }
            //                                // Spacing 2 is bad
            //                                else if (!spacingMassRange[2].Contains(spacing2))
            //                                {
            //                                    pair.peaks[j + 3, k] = null;
            //                                    peakCount--;
            //                                }
            //                            }
            //                            // Peaks 2,3 & 5
            //                            else if (peak5 != null)
            //                            {
            //                                peak3NeutralMass = Mass.MassFromMz(peak5.X, charge);
            //                                spacing2 = peak3NeutralMass - peak2NeutralMass;
            //                                double extendedTheoSpacing = spacingMassRange[2].Mean + spacingMassRange[3].Mean;
            //                                MassRange doubleSpacing = new MassRange(extendedTheoSpacing, new MassTolerance(MassToleranceType.DA, extendedTheoSpacing * SPACINGPERCENTAGEERROR));

            //                                // Spacing 1 is bad
            //                                if (!spacingMassRange[1].Contains(spacing1))
            //                                {
            //                                    pair.peaks[j + 1, k] = null;
            //                                    peakCount--;

            //                                    // Spacings 1 & 2 are bad
            //                                    if (!doubleSpacing.Contains(spacing2))
            //                                    {
            //                                        pair.peaks[j + 2, k] = null;
            //                                        peakCount--;
            //                                        pair.peaks[j + 4, k] = null;
            //                                        peakCount--;
            //                                    }
            //                                }
            //                                // Spacing 2 is bad
            //                                else if (!doubleSpacing.Contains(spacing2))
            //                                {
            //                                    pair.peaks[j + 4, k] = null;
            //                                    peakCount--;
            //                                }
            //                            }
            //                            // Peaks 2,3 & 6
            //                            else if (peak6 != null)
            //                            {
            //                                peak3NeutralMass = Mass.MassFromMz(peak6.X, charge);
            //                                spacing2 = peak3NeutralMass - peak2NeutralMass;
            //                                double extendedTheoSpacing = spacingMassRange[2].Mean + spacingMassRange[3].Mean + spacingMassRange[4].Mean;
            //                                MassRange tripleSpacing = new MassRange(extendedTheoSpacing, new MassTolerance(MassToleranceType.DA, extendedTheoSpacing * SPACINGPERCENTAGEERROR));

            //                                // Spacing 1 is bad
            //                                if (!spacingMassRange[1].Contains(spacing1))
            //                                {
            //                                    pair.peaks[j + 1, k] = null;
            //                                    peakCount--;

            //                                    // Spacings 1 & 2 are bad
            //                                    if (!tripleSpacing.Contains(spacing2))
            //                                    {
            //                                        pair.peaks[j + 2, k] = null;
            //                                        peakCount--;
            //                                        pair.peaks[j + 5, k] = null;
            //                                        peakCount--;
            //                                    }
            //                                }
            //                                // Spacing 2 is bad
            //                                else if (!tripleSpacing.Contains(spacing2))
            //                                {
            //                                    pair.peaks[j + 5, k] = null;
            //                                    peakCount--;
            //                                }
            //                            }
            //                        }
            //                        else if (peak4 != null)
            //                        {
            //                            peak2NeutralMass = Mass.MassFromMz(peak4.X, charge);
            //                            spacing1 = peak2NeutralMass - peak1NeutralMass;
            //                            double extendedTheoSpacing = spacingMassRange[1].Mean + spacingMassRange[2].Mean;
            //                            MassRange doubleSpacing = new MassRange(extendedTheoSpacing, new MassTolerance(MassToleranceType.DA, extendedTheoSpacing * SPACINGPERCENTAGEERROR));

            //                            // Peaks 2,4 & 5
            //                            if (peak5 != null)
            //                            {
            //                                peak3NeutralMass = Mass.MassFromMz(peak5.X, charge);
            //                                spacing2 = peak3NeutralMass - peak2NeutralMass;

            //                                // Spacing 1 is bad
            //                                if (!doubleSpacing.Contains(spacing1))
            //                                {
            //                                    pair.peaks[j + 1, k] = null;
            //                                    peakCount--;

            //                                    // Spacings 1 & 2 are bad
            //                                    if (!spacingMassRange[3].Contains(spacing2))
            //                                    {
            //                                        pair.peaks[j + 3, k] = null;
            //                                        peakCount--;
            //                                        pair.peaks[j + 4, k] = null;
            //                                        peakCount--;
            //                                    }
            //                                }
            //                                // Spacing 2 is bad
            //                                else if (!spacingMassRange[3].Contains(spacing2))
            //                                {
            //                                    pair.peaks[j + 4, k] = null;
            //                                    peakCount--;
            //                                }
            //                            }
            //                            // Peaks 2,4 & 6
            //                            else if (peak6 != null)
            //                            {
            //                                peak3NeutralMass = Mass.MassFromMz(peak6.X, charge);
            //                                spacing2 = peak3NeutralMass - peak2NeutralMass;
            //                                double extendedTheoSpacing1 = spacingMassRange[3].Mean + spacingMassRange[4].Mean;
            //                                MassRange doubleSpacing1 = new MassRange(extendedTheoSpacing1, new MassTolerance(MassToleranceType.DA, extendedTheoSpacing1 * SPACINGPERCENTAGEERROR));

            //                                // Spacing 1 is bad
            //                                if (!doubleSpacing.Contains(spacing1))
            //                                {
            //                                    pair.peaks[j + 1, k] = null;
            //                                    peakCount--;

            //                                    // Spacings 1 & 2 are bad
            //                                    if (!doubleSpacing1.Contains(spacing2))
            //                                    {
            //                                        pair.peaks[j + 3, k] = null;
            //                                        peakCount--;
            //                                        pair.peaks[j + 5, k] = null;
            //                                        peakCount--;
            //                                    }
            //                                }
            //                                // Spacing 2 is bad
            //                                else if (!doubleSpacing1.Contains(spacing2))
            //                                {
            //                                    pair.peaks[j + 5, k] = null;
            //                                    peakCount--;
            //                                }
            //                            }
            //                        }
            //                        // Peaks 2,5 & 6
            //                        else if (peak5 != null)
            //                        {
            //                            peak2NeutralMass = Mass.MassFromMz(peak5.X, charge);
            //                            peak3NeutralMass = Mass.MassFromMz(peak6.X, charge);
            //                            spacing1 = peak2NeutralMass - peak1NeutralMass;
            //                            spacing2 = peak3NeutralMass - peak2NeutralMass;
            //                            double extendedTheoSpacing = spacingMassRange[1].Mean + spacingMassRange[2].Mean + spacingMassRange[3].Mean;
            //                            MassRange tripleSpacing = new MassRange(extendedTheoSpacing, new MassTolerance(MassToleranceType.DA, extendedTheoSpacing * SPACINGPERCENTAGEERROR));

            //                            // Spacing 1 is bad
            //                            if (!tripleSpacing.Contains(spacing1))
            //                            {
            //                                pair.peaks[j + 1, k] = null;
            //                                peakCount--;

            //                                // Spacings 1 & 2 are bad
            //                                if (!spacingMassRange[4].Contains(spacing2))
            //                                {
            //                                    pair.peaks[j + 4, k] = null;
            //                                    peakCount--;
            //                                    pair.peaks[j + 5, k] = null;
            //                                    peakCount--;
            //                                }
            //                            }
            //                            // Spacing 2 is bad
            //                            else if (!spacingMassRange[4].Contains(spacing2))
            //                            {
            //                                pair.peaks[j + 5, k] = null;
            //                                peakCount--;
            //                            }
            //                        }
            //                    }
            //                    else if (peak3 != null)
            //                    {
            //                        peak1NeutralMass = Mass.MassFromMz(peak3.X, charge);

            //                        if (peak4 != null)
            //                        {
            //                            peak2NeutralMass = Mass.MassFromMz(peak4.X, charge);
            //                            spacing1 = peak2NeutralMass - peak1NeutralMass;

            //                            // Peaks 3,4 & 5
            //                            if (peak5 != null)
            //                            {
            //                                peak3NeutralMass = Mass.MassFromMz(peak5.X, charge);
            //                                spacing2 = peak3NeutralMass - peak2NeutralMass;

            //                                // Spacing 1 is bad
            //                                if (!spacingMassRange[2].Contains(spacing1))
            //                                {
            //                                    pair.peaks[j + 2, k] = null;
            //                                    peakCount--;

            //                                    // Spacings 1 & 2 are bad
            //                                    if (!spacingMassRange[3].Contains(spacing2))
            //                                    {
            //                                        pair.peaks[j + 3, k] = null;
            //                                        peakCount--;
            //                                        pair.peaks[j + 4, k] = null;
            //                                        peakCount--;
            //                                    }
            //                                }
            //                                // Spacing 2 is bad
            //                                else if (!spacingMassRange[3].Contains(spacing2))
            //                                {
            //                                    pair.peaks[j + 4, k] = null;
            //                                    peakCount--;
            //                                }
            //                            }
            //                            // Peaks 3,4 & 6
            //                            else if (peak6 != null)
            //                            {
            //                                peak3NeutralMass = Mass.MassFromMz(peak6.X, charge);
            //                                spacing2 = peak3NeutralMass - peak2NeutralMass;
            //                                double extendedTheoSpacing = spacingMassRange[3].Mean + spacingMassRange[4].Mean;
            //                                MassRange doubleSpacing = new MassRange(extendedTheoSpacing, new MassTolerance(MassToleranceType.DA, extendedTheoSpacing * SPACINGPERCENTAGEERROR));

            //                                // Spacing 1 is bad
            //                                if (!spacingMassRange[2].Contains(spacing1))
            //                                {
            //                                    pair.peaks[j + 2, k] = null;
            //                                    peakCount--;

            //                                    // Spacings 1 & 2 are bad
            //                                    if (!doubleSpacing.Contains(spacing2))
            //                                    {
            //                                        pair.peaks[j + 3, k] = null;
            //                                        peakCount--;
            //                                        pair.peaks[j + 5, k] = null;
            //                                        peakCount--;
            //                                    }
            //                                }
            //                                // Spacing 2 is bad
            //                                else if (!doubleSpacing.Contains(spacing2))
            //                                {
            //                                    pair.peaks[j + 5, k] = null;
            //                                    peakCount--;
            //                                }
            //                            }
            //                        }
            //                        // Peaks 3,5 & 6
            //                        else if (peak5 != null)
            //                        {
            //                            peak2NeutralMass = Mass.MassFromMz(peak5.X, charge);
            //                            peak3NeutralMass = Mass.MassFromMz(peak6.X, charge);
            //                            spacing1 = peak2NeutralMass - peak1NeutralMass;
            //                            spacing2 = peak3NeutralMass - peak2NeutralMass;
            //                            double extendedTheoSpacing = spacingMassRange[2].Mean + spacingMassRange[3].Mean;
            //                            MassRange doubleSpacing = new MassRange(extendedTheoSpacing, new MassTolerance(MassToleranceType.DA, extendedTheoSpacing * SPACINGPERCENTAGEERROR));

            //                            // Spacing 1 is bad
            //                            if (!doubleSpacing.Contains(spacing1))
            //                            {
            //                                pair.peaks[j + 2, k] = null;
            //                                peakCount--;

            //                                // Spacings 1 & 2 are bad
            //                                if (!spacingMassRange[4].Contains(spacing2))
            //                                {
            //                                    pair.peaks[j + 4, k] = null;
            //                                    peakCount--;
            //                                    pair.peaks[j + 5, k] = null;
            //                                    peakCount--;
            //                                }
            //                            }
            //                            // Spacing 2 is bad
            //                            else if (!spacingMassRange[4].Contains(spacing2))
            //                            {
            //                                pair.peaks[j + 5, k] = null;
            //                                peakCount--;
            //                            }
            //                        }
            //                    }
            //                    // Peaks 4,5 & 6
            //                    else if (peak4 != null)
            //                    {
            //                        peak1NeutralMass = Mass.MassFromMz(peak4.X, charge);
            //                        peak2NeutralMass = Mass.MassFromMz(peak5.X, charge);
            //                        peak3NeutralMass = Mass.MassFromMz(peak6.X, charge);
            //                        spacing1 = peak2NeutralMass - peak1NeutralMass;
            //                        spacing2 = peak3NeutralMass - peak2NeutralMass;

            //                        // Spacing 1 is bad
            //                        if (!spacingMassRange[3].Contains(spacing1))
            //                        {
            //                            pair.peaks[j + 3, k] = null;
            //                            peakCount--;

            //                            // Spacings 1 & 2 are bad
            //                            if (!spacingMassRange[4].Contains(spacing2))
            //                            {
            //                                pair.peaks[j + 4, k] = null;
            //                                peakCount--;
            //                                pair.peaks[j + 5, k] = null;
            //                                peakCount--;
            //                            }
            //                        }
            //                        // Spacing 2 is bad
            //                        else if (!spacingMassRange[4].Contains(spacing2))
            //                        {
            //                            pair.peaks[j + 5, k] = null;
            //                            peakCount--;
            //                        }
            //                    }
            //                }
            //                else if (Form1.NOISEBANDCAP && peakCount < 3)
            //                {
            //                    pair.peaks[j, k] = null;
            //                    pair.peaks[j + 1, k] = null;
            //                    pair.peaks[j + 2, k] = null;
            //                    pair.peaks[j + 3, k] = null;
            //                    pair.peaks[j + 4, k] = null;
            //                    pair.peaks[j + 5, k] = null;
            //                }

            //                if (peakCount < peaksNeeded)
            //                {
            //                    pair.peaks[j, k] = null;
            //                    pair.peaks[j + 1, k] = null;
            //                    pair.peaks[j + 2, k] = null;
            //                    pair.peaks[j + 3, k] = null;
            //                    pair.peaks[j + 4, k] = null;
            //                    pair.peaks[j + 5, k] = null;
            //                }
            //            }
            //        }
            //    }
            //}

            List<Pair> spacingsCheckedList = new List<Pair>();
            List<Pair> completePairList = null;
            completePairs.TryGetValue(rawFile, out completePairList);

            if (numIsotopologues > 1)
            {
                foreach (Pair currentPair in pairs)
                {
                    for (int c = 0; c < numClusters; c++)
                    {
                        for (int i = 0; i < numIsotopes; i++)
                        {
                            if (currentPair.complete[c, i])
                            {
                                completePairList = addPairToList(completePairList, currentPair, c, i, rawFile);
                            }
                            else if (currentPair.peakCount[c, i] >= peaksNeeded)
                            {
                                spacingsCheckedList = addPairToList(spacingsCheckedList, currentPair, c, i, rawFile);
                            }
                        }
                    }
                }
                completePairs[rawFile] = completePairList;
                allPairs[rawFile] = spacingsCheckedList;
            }
            else
            {
                foreach (Pair currentPair in pairs)
                {
                    for (int i = 0; i < numIsotopes; i++)
                    {
                        if (currentPair.complete[0, i])
                        {
                            completePairList = addPairToList(completePairList, currentPair, 0, i, rawFile);
                        }
                        else if (currentPair.peakCount[0, i] >= peaksNeeded)
                        {
                            spacingsCheckedList = addPairToList(completePairList, currentPair, 0, i, rawFile);
                        }
                    }
                }
                completePairs[rawFile] = completePairList;
                allPairs[rawFile] = spacingsCheckedList;
            }
            
            // When noise-band capping not enabled, add pairs to complete list
            //if (!Form1.NOISEBANDCAP)
            //{
            //    spacingsCheckedList = new List<Pair>();
            //    foreach (Pair currentPair in pairs)
            //    {
            //        if (currentPair.totalPeakCount > 0)
            //        {
            //            spacingsCheckedList.Add(currentPair);
            //        }
            //    }
            //    allPairs[rawFile] = spacingsCheckedList;
            //    completePairs[rawFile] = spacingsCheckedList;
            //}
            //else
            //{
            //    completePairs.TryGetValue(rawFile, out completePairList);
            //    spacingsCheckedList = new List<Pair>();
            //    if (numIsotopologues > 1)
            //    {
            //        foreach (Pair currentPair in pairs)
            //        {
            //            if (currentPair.totalPeakCount > 0)
            //            {
            //                spacingsCheckedList.Add(currentPair);
            //            }
                        
            //            for (int c = 0; c < numClusters; c++)
            //            {
            //                for (int i = 0; i < numIsotopes; i++)
            //                {
            //                    if (currentPair.complete[c, i])
            //                    {
            //                        addToCompleteList(currentPair, c, i, rawFile);
            //                    }
            //                }
            //            }
            //        }
            //        allPairs[rawFile] = spacingsCheckedList;
            //    }
            //    else
            //    {
            //        foreach (Pair currentPair in pairs)
            //        {
            //            if (currentPair.totalPeakCount > 0)
            //            {
            //                spacingsCheckedList.Add(currentPair);
            //            }
            //            for (int i = 0; i < numIsotopes; i++)
            //            {
            //                if (currentPair.complete[0, i])
            //                {
            //                    addToCompleteList(currentPair, 0, i, rawFile);
            //                }
            //            }
            //        }
            //        allPairs[rawFile] = spacingsCheckedList;
            //    }                
            //}
        }

        public void checkIsotopeDistribution(MSDataFile rawFile, List<Pair> pairs)
        {            
            List<Pair> goodIsotopeDistributionPairs = new List<Pair>();

            if (numIsotopologues > 1)
            {
                foreach (Pair pair in pairs)
                {
                    for (int c = 0; c < numClusters; c++)
                    {
                        if (numIsotopes > 1 && pair.checkIsotopeDistribution(c))
                        {
                            for (int j = 0; j < numIsotopes; j++)
                            {
                                if (pair.peakCount[c, j] < peaksNeeded)
                                {
                                    int channelIndex = c * numIsotopologues;
                                    for (int i = channelIndex; i < channelIndex + numIsotopologues; i++)
                                    {
                                        pair.peaks[i, j] = null;
                                    }
                                }
                            }
                        }
                    }
                    if (pair.totalPeakCount > 0)
                    {
                        goodIsotopeDistributionPairs.Add(pair);
                    }
                }
            }
            else
            {
                foreach (Pair pair in pairs)
                {
                    for (int c = 0; c < numChannels; c++)
                    {
                        pair.checkIsotopeDistribution(c);
                    }

                    for (int j = 0; j < numIsotopes; j++)
                    {
                        if (pair.peakCount[0, j] < peaksNeeded)
                        {
                            for (int c = 0; c < numChannels; c++)
                            {
                                pair.peaks[c, j] = null;
                            }
                        }
                    }

                    if (pair.totalPeakCount > 0)
                    {
                        goodIsotopeDistributionPairs.Add(pair);
                    }
                }
            }
            allPairs[rawFile] = goodIsotopeDistributionPairs;
        }

        public void sortPairs(MSDataFile rawFile)
        {
            PeakPattern[,] pairPatterns;
            List<Pair> associatedPairs = null;

            if (allPairs.TryGetValue(rawFile, out associatedPairs))
            {
                // NeuCode quantification
                if (numIsotopologues > 1)
                {
                    for (int c = 0; c < numClusters; c++)
                    {
                        int channelIndex = c * numIsotopologues;
                        List<Pair> newAllPairs = new List<Pair>();

                        // Only look to sort pairs if there are 3 or more incomplete sets
                        if ((countAllIsotopes[c] - countCompleteIsotopes[c]) >= 3)
                        {
                            int numPeaks = numIsotopologues - peaksNeeded;
                            int patternPossibilities = numIsotopologues * numPeaks;

                            pairPatterns = new PeakPattern[numIsotopologues, patternPossibilities];

                            for (int i = 0; i < numIsotopologues; i++)
                            {
                                for (int j = 0; j < patternPossibilities; j++)
                                {
                                    pairPatterns[i, j] = new PeakPattern(i, j, this);
                                }
                            }

                            foreach (Pair pair in associatedPairs)
                            {
                                for (int j = 0; j < numIsotopes; j++)
                                {
                                    if (pair.peakCount[c, j] >= peaksNeeded && pair.peakCount[c, j] < numIsotopologues && pair.peakPattern[c, j] >= 0)
                                    {
                                        Pair cleanedPair = new Pair(pair.parent, pair.rawFile, pair.scanNumber, pair.injectionTime, pair.retentionTime);
                                        for (int i = 0; i < numChannels; i++)
                                        {
                                            for (int k = 0; k < numIsotopes; k++)
                                            {
                                                if (k == j && i >= channelIndex && i < channelIndex + numIsotopologues)
                                                {
                                                    cleanedPair.peaks[i, k] = pair.peaks[i, k];
                                                }
                                            }
                                        }
                                        pairPatterns[pair.peakCount[c, j], pair.peakPattern[c, j]].PatternPairs.Add(cleanedPair);
                                    }
                                }
                            }

                            // Find maximum pattern score and perform normalization
                            double maximumScore = 0;
                            for (int i = 0; i < numIsotopologues; i++)
                            {
                                for (int j = 0; j < patternPossibilities; j++)
                                {
                                    if (pairPatterns[i, j].RawScore > maximumScore)
                                    {
                                        maximumScore = pairPatterns[i, j].RawScore;
                                    }
                                }
                            }

                            maximumScore += 1;

                            for (int i = 0; i < numIsotopologues; i++)
                            {
                                for (int j = 0; j < patternPossibilities; j++)
                                {
                                    pairPatterns[i, j].setNormalizedScore(maximumScore);
                                }
                            }

                            // Find best pattern(s)

                            double bestScore = 0;
                            double secondBestScore = 0;
                            int bestNumPeaks = 0;
                            int bestPatternPossibilities = 0;

                            for (int i = 0; i < numIsotopologues; i++)
                            {
                                for (int j = 0; j < patternPossibilities; j++)
                                {
                                    if (pairPatterns[i, j].NormalizedScore > bestScore)
                                    {
                                        bestNumPeaks = i;
                                        bestPatternPossibilities = j;
                                        secondBestScore = bestScore;
                                        bestScore = pairPatterns[i, j].NormalizedScore;
                                    }
                                    else if (pairPatterns[i, j].NormalizedScore > secondBestScore)
                                    {
                                        secondBestScore = pairPatterns[i, j].NormalizedScore;
                                    }
                                }
                            }

                            // If one pattern has 3 or more measurements

                            if (bestNumPeaks > 0 && bestScore >= 3 && (bestScore / secondBestScore) > 3.0)
                            {
                                missingChannelPattern[c, 0] = bestNumPeaks;
                                missingChannelPattern[c, 1] = bestPatternPossibilities;

                                foreach (Pair pair in pairPatterns[bestNumPeaks, bestPatternPossibilities].PatternPairs)
                                {
                                    newAllPairs.Add(pair);
                                }
                            }

                            // If the best pattern has fewer than 3 measurements, group similar patterns together

                            //else
                            //{
                            //    if (bestNumPeaks == 2)
                            //    {
                            //        int altNumPeaks = 3;
                            //        int altPattern1 = -1;
                            //        int altPattern2 = -1;

                            //        switch (bestPatternPossibilities)
                            //        {
                            //            case 0:
                            //                altPattern1 = 0;
                            //                altPattern2 = 1;
                            //                break;
                            //            case 1:
                            //                altPattern1 = 0;
                            //                altPattern2 = 2;
                            //                break;
                            //            case 2:
                            //                altPattern1 = 1;
                            //                altPattern2 = 2;
                            //                break;
                            //            case 3:
                            //                altPattern1 = 0;
                            //                altPattern2 = 3;
                            //                break;
                            //            case 4:
                            //                altPattern1 = 1;
                            //                altPattern2 = 3;
                            //                break;
                            //            case 5:
                            //                altPattern1 = 2;
                            //                altPattern2 = 3;
                            //                break;
                            //        }

                            //        foreach (Pair pair in pairPatterns[bestNumPeaks, bestPatternPossibilities].PatternPairs)
                            //        {
                            //            newAllPairs.Add(pair);
                            //        }

                            //        foreach (Pair pair in pairPatterns[altNumPeaks, altPattern1].PatternPairs)
                            //        {
                            //            newAllPairs.Add(pair);
                            //        }

                            //        foreach (Pair pair in pairPatterns[altNumPeaks, altPattern2].PatternPairs)
                            //        {
                            //            newAllPairs.Add(pair);
                            //        }
                            //    }
                            //    else if (bestNumPeaks == 3)
                            //    {
                            //        int altNumPeaks = 2;
                            //        int altPattern1 = -1;
                            //        int altPattern2 = -1;
                            //        int altPattern3 = -1;

                            //        switch (bestPatternPossibilities)
                            //        {
                            //            case 0:
                            //                altPattern1 = 0;
                            //                altPattern2 = 1;
                            //                altPattern3 = 3;
                            //                break;
                            //            case 1:
                            //                altPattern1 = 0;
                            //                altPattern2 = 2;
                            //                altPattern3 = 4;
                            //                break;
                            //            case 2:
                            //                altPattern1 = 1;
                            //                altPattern2 = 2;
                            //                altPattern3 = 5;
                            //                break;
                            //            case 3:
                            //                altPattern1 = 3;
                            //                altPattern2 = 4;
                            //                altPattern3 = 5;
                            //                break;
                            //        }

                            //        foreach (Pair pair in pairPatterns[bestNumPeaks, bestPatternPossibilities].PatternPairs)
                            //        {
                            //            newAllPairs.Add(pair);
                            //        }

                            //        foreach (Pair pair in pairPatterns[altNumPeaks, altPattern1].PatternPairs)
                            //        {
                            //            newAllPairs.Add(pair);
                            //        }

                            //        foreach (Pair pair in pairPatterns[altNumPeaks, altPattern2].PatternPairs)
                            //        {
                            //            newAllPairs.Add(pair);
                            //        }

                            //        foreach (Pair pair in pairPatterns[altNumPeaks, altPattern3].PatternPairs)
                            //        {
                            //            newAllPairs.Add(pair);
                            //        }
                            //    }
                            //}

                            if (newAllPairs != null && newAllPairs.Count >= 3)
                            {
                                allPairs[rawFile] = newAllPairs;
                            }
                        }
                    }
                }
                // Non-NeuCode quantification
                else
                {
                    List<Pair> newAllPairs = new List<Pair>();

                    // Only look to sort pairs if there are 3 or more incomplete sets
                    if ((countAllIsotopes[0] - countCompleteIsotopes[0]) >= 3)
                    {
                        int numPeaks = numChannels - peaksNeeded;
                        int patternPossibilities = numChannels * numPeaks;

                        pairPatterns = new PeakPattern[numChannels, patternPossibilities];

                        for (int i = 0; i < numChannels; i++)
                        {
                            for (int j = 0; j < patternPossibilities; j++)
                            {
                                pairPatterns[i, j] = new PeakPattern(i, j, this);
                            }
                        }

                        foreach (Pair pair in associatedPairs)
                        {
                            for (int j = 0; j < numIsotopes; j++)
                            {
                                if (pair.peakCount[0, j] >= peaksNeeded && pair.peakCount[0, j] < numChannels && pair.peakPattern[0, j] >= 0)
                                {
                                    Pair cleanedPair = new Pair(pair.parent, pair.rawFile, pair.scanNumber, pair.injectionTime, pair.retentionTime);
                                    for (int i = 0; i < numChannels; i++)
                                    {
                                        for (int k = 0; k < numIsotopes; k++)
                                        {
                                            if (k == j && i >= 0 && i < numChannels)
                                            {
                                                cleanedPair.peaks[i, k] = pair.peaks[i, k];
                                            }
                                        }
                                    }
                                    pairPatterns[pair.peakCount[0, j], pair.peakPattern[0, j]].PatternPairs.Add(cleanedPair);
                                }
                            }
                        }

                        // Find maximum pattern score and perform normalization
                        double maximumScore = 0;
                        for (int i = 0; i < numChannels; i++)
                        {
                            for (int j = 0; j < patternPossibilities; j++)
                            {
                                if (pairPatterns[i, j].RawScore > maximumScore)
                                {
                                    maximumScore = pairPatterns[i, j].RawScore;
                                }
                            }
                        }

                        maximumScore += 1;

                        for (int i = 0; i < numChannels; i++)
                        {
                            for (int j = 0; j < patternPossibilities; j++)
                            {
                                pairPatterns[i, j].setNormalizedScore(maximumScore);
                            }
                        }

                        // Find best pattern(s)

                        double bestScore = 0;
                        int bestNumPeaks = 0;
                        int bestPatternPossibilities = 0;

                        for (int i = 0; i < numChannels; i++)
                        {
                            for (int j = 0; j < patternPossibilities; j++)
                            {
                                if (pairPatterns[i, j].NormalizedScore > bestScore)
                                {
                                    bestNumPeaks = i;
                                    bestPatternPossibilities = j;
                                    bestScore = pairPatterns[i, j].NormalizedScore;
                                }
                            }
                        }

                        // If one pattern has 3 or more measurements

                        if (bestNumPeaks > 0 && bestScore >= 3)
                        {
                            missingChannelPattern[0, 0] = bestNumPeaks;
                            missingChannelPattern[0, 1] = bestPatternPossibilities;

                            foreach (Pair pair in pairPatterns[bestNumPeaks, bestPatternPossibilities].PatternPairs)
                            {
                                newAllPairs.Add(pair);
                            }
                        }

                        if (newAllPairs != null && newAllPairs.Count >= 3)
                        {
                            allPairs[rawFile] = newAllPairs;
                        }
                    }
                }
            }
        }

        public List<Pair> addPairToList(List<Pair> initialList, Pair pair, int cluster, int isotope, MSDataFile rawFile)
        {
            List<Pair> updatedList = new List<Pair>();
            if (initialList == null) initialList = new List<Pair>();
            updatedList = initialList;
            Pair pairToAdd = new Pair(this, rawFile, pair.scanNumber, pair.injectionTime, pair.retentionTime);
            ILabeledPeak[,] peaksToAdd = new ILabeledPeak[numChannels, numIsotopes];

            if (numIsotopologues > 1)
            {
                int channelIndex = cluster * numIsotopologues;
                for (int i = channelIndex; i < channelIndex + numIsotopologues; i++)
                {
                    peaksToAdd[i, isotope] = pair.peaks[i, isotope];
                }

                if (updatedList.Count == 0)
                {
                    pairToAdd.peaks = peaksToAdd;
                    updatedList.Add(pairToAdd);
                }
                else
                {
                    bool found = false;
                    foreach (Pair oldPair in updatedList)
                    {
                        if (oldPair.scanNumber == pairToAdd.scanNumber)
                        {
                            for (int i = channelIndex; i < channelIndex + numIsotopologues; i++)
                            {
                                oldPair.peaks[i, isotope] = peaksToAdd[i, isotope];
                            }
                            found = true;
                        }
                    }
                    if (!found)
                    {
                        pairToAdd.peaks = peaksToAdd;
                        updatedList.Add(pairToAdd);
                    }
                }
            }
            else
            {
                for (int i = 0; i < numChannels; i++)
                {
                    peaksToAdd[i, isotope] = pair.peaks[i, isotope];
                }

                if (updatedList.Count == 0)
                {
                    pairToAdd.peaks = peaksToAdd;
                    updatedList.Add(pairToAdd);
                }
                else
                {
                    bool found = false;
                    foreach (Pair oldPair in updatedList)
                    {
                        if (oldPair.scanNumber == pairToAdd.scanNumber)
                        {
                            for (int i = 0; i < numChannels; i++)
                            {
                                oldPair.peaks[i, isotope] = peaksToAdd[i, isotope];
                            }
                            found = true;
                        }
                    }
                    if (!found)
                    {
                        pairToAdd.peaks = peaksToAdd;
                        updatedList.Add(pairToAdd);
                    }
                }
            }
            return updatedList;
        }

        /* To missing channels of pairs that are not coalesced, an noise-based intensity is added to permit quantification
         */
        public void applyNoise(MSDataFile rawFile)
        {
            //if (sequence.Equals("TVFGLSDSMSK"))
            //{
            //    int found = 0;
            //}
            
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

            if (numIsotopologues > 1)
            {
                foreach (Pair pair in all)
                {
                    for (int j = 0; j < numIsotopes; j++)
                    {
                        for (int c = 0; c < numClusters; c++)
                        {
                            channelIndex = c * numIsotopologues;
                            //// For complete pairs
                            //if (pair.complete[c, j])
                            //{
                            //    // If the pair falls above the peptide's coalescence threshold, set peaks to null
                            //    if (coalescenceDetected && pair.maxIntensity[c, j] > coalescenceIntensity)
                            //    {
                            //        for (int m = channelIndex; m < channelIndex + numIsotopologues; m++)
                            //        {
                            //            pair.peaks[m, j] = null;
                            //        }
                            //    }
                            //    // Otherwise, add pair to complete list
                            //    else
                            //    {
                            //        // First check to see if the peptide is already on the complete list
                            //        bool found = false;
                            //        int SN = pair.scanNumber;

                            //        foreach (Pair noNBC in completeOnly)
                            //        {
                            //            if (noNBC.scanNumber == SN)
                            //            {
                            //                // If found, update peaks for that pair
                            //                for (int m = channelIndex; m < channelIndex + numIsotopologues; m++)
                            //                {
                            //                    noNBC.peaks[m, j] = pair.peaks[m, j];
                            //                }
                            //                found = true;
                            //            }
                            //        }

                            //        // If not found, add pair
                            //        if (!found)
                            //        {
                            //            Pair noNBCPair = new Pair(this, rawFile, pair.scanNumber, pair.injectionTime, pair.retentionTime);
                            //            for (int m = channelIndex; m < channelIndex + numIsotopologues; m++)
                            //            {
                            //                noNBCPair.peaks[m, j] = pair.peaks[m, j];
                            //            }
                            //            completeOnly.Add(noNBCPair);
                            //        }
                            //    }
                            //}
                            // For incomplete pairs
                            if (!pair.complete[c,j] && pair.peakCount[c, j] >= peaksNeeded)
                            {
                                // Check for coalescence
                                if (coalescenceDetected && pair.maxClusterIntensity[c, j] > coalescenceIntensity)
                                {
                                    if (pair.maxClusterIntensity[c, j] > coalescenceIntensity)
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
            else
            {
                foreach (Pair pair in all)
                {
                    for (int j = 0; j < numIsotopes; j++)
                    {
                        // For complete pairs
                        if (pair.complete[0, j])
                        {
                            // First check to see if the peptide is already on the complete list
                            bool found = false;
                            int SN = pair.scanNumber;

                            foreach (Pair noNBC in completeOnly)
                            {
                                if (noNBC.scanNumber == SN)
                                {
                                    // If found, update peaks for that pair
                                    for (int m = 0; m < numChannels; m++)
                                    {
                                        noNBC.peaks[m, j] = pair.peaks[m, j];
                                    }
                                    found = true;
                                }
                            }

                            // If not found, add pair
                            if (!found)
                            {
                                Pair noNBCPair = new Pair(this, rawFile, pair.scanNumber, pair.injectionTime, pair.retentionTime);
                                for (int m = 0; m < numChannels; m++)
                                {
                                    noNBCPair.peaks[m, j] = pair.peaks[m, j];
                                }
                                completeOnly.Add(noNBCPair);
                            }
                        }

                        // For incomplete pairs
                        else if (pair.peakCount[0, j] >= peaksNeeded)
                        {
                            // Apply noise to missing channels
                            for (int m = 0; m < numChannels; m++)
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
                //if (Form1.NEUCODE_SIXPLEX_ARG && Form1.CORRECTARGPROLINE && countResidues('R', sequence) > 0 && countResidues('P', sequence) > 0)
                //{
                //    int correction = correctIsotopeDistributions(lowerResolutionScan, rawFile);
                //    theoMasses[4, 0] += correction * (5 * (Constants.CARBON13 - Constants.CARBON));
                //    theoMasses[5, 0] += correction * (5 * (Constants.CARBON13 - Constants.CARBON));
                //}
                //if (Form1.NEUCODE_SIXPLEX_LEU && Form1.CORRECTLEUDLOSS && countResidues('L', sequence) > 0)
                //{
                //    int correction = correctIsotopeDistributions(lowerResolutionScan, rawFile);
                //    theoMasses[4, 0] -= correction * (Constants.DEUTERIUM - Constants.HYDROGEN);
                //    theoMasses[5, 0] -= correction * (Constants.DEUTERIUM - Constants.HYDROGEN);
                //}
                //if (Form1.NEUCODE_DUPLEX_LEU7_18MDA && Form1.CORRECTLEUNLOSS && countResidues('L', sequence) > 0)
                //{
                //    int correction = correctIsotopeDistributions(lowerResolutionScan, rawFile);
                //    theoMasses[0, 0] -= correction * (Constants.NITROGEN15 - Constants.NITROGEN);
                //}
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
                ILabeledPeak[] peaks = new ILabeledPeak[numChannels];
                int peaksCount = 0;
                
                //if (Form1.SILAC_DUPLEX_LEUCN && Form1.CORRECTLEUNLOSS)
                //{
                //    int correction = correctIsotopeDistributions(precursorScan, rawFile);
                //    theoMasses[1, 0] -= correction * (Constants.NITROGEN15 - Constants.NITROGEN);
                //}
                for (int i = 0; i < numChannels; i++)
                {
                    peaks[i] = largestPeak(theoMasses[i, 0], precursorScan, firstSearchMassRange, rawFile);
                    if (peaks[i] != null) peaksCount++;
                }
                if (peaksCount == numChannels)
                {
                    for (int i = 0; i < numChannels; i++)
                    {
                        precursorPPM = MassTolerance.GetTolerance(Mass.MassFromMz(peaks[i].X, charge), theoMasses[i, 0], MassToleranceType.PPM);
                        ppm = new PrecursorPPM(charge, this.sequence, best.EValue, precursorPPM);
                        ppms.Add(ppm);
                    }    
                }
            }             
        }

        //public int correctIsotopeDistributions(MSDataScan current, MSDataFile rawFile)
        //{
        //    int numLabels;
        //    if (Form1.SILAC_DUPLEX_LEUCN || Form1.NEUCODE_DUPLEX_LEU7_18MDA || Form1.NEUCODE_SIXPLEX_LEU)
        //    {
        //        numLabels = countResidues('L', sequence) + 1;
        //    }
        //    else if (Form1.NEUCODE_SIXPLEX_ARG && countResidues('R', sequence) > 0)
        //    {
        //        numLabels = countResidues('P', sequence) + 1;
        //    }
        //    else
        //    {
        //        numLabels = 0;
        //    }

        //    double[,] correctedMasses = new double[numLabels, numIsotopes];// Rows shift theoretical masses to check for isotope conversion, columns are isotopes
        //    double[,] correctedIntensities = new double[numLabels, numIsotopes]; // Rows shift theoretical masses to check for isotope conversion, columns are isotopes
        //    int bestMonoMass = 0; // 0 = uncorrected, 1 = corrected for 1 label loss, etc.
        //    double uncorrectedMonoMass;
        //    if ((Form1.SILAC_DUPLEX_LEUCN || Form1.NEUCODE_DUPLEX_LEU7_18MDA) && Form1.CORRECTLEUNLOSS && countResidues('L', sequence) > 0)
        //    {
        //        // Fill first row
        //        if (Form1.SILAC_DUPLEX_LEUCN)
        //        {
        //            uncorrectedMonoMass = theoMasses[1, 0];
        //        }
        //        else if (Form1.NEUCODE_DUPLEX_LEU7_18MDA)
        //        {
        //            uncorrectedMonoMass = (theoMasses[0, 0] + theoMasses[1, 0]) / 2.0;
        //        }
        //        else
        //        {
        //            uncorrectedMonoMass = theoMasses[0, 0];
        //        }

        //        // Fill rest of the array
        //        double mass;
        //        for (int i = 0; i < numLabels; i++)
        //        {
        //            mass = uncorrectedMonoMass - (i * (Constants.NITROGEN15 - Constants.NITROGEN));
        //            for (int j = 0; j < numIsotopes; j++)
        //            {
        //                correctedMasses[i, j] = mass + (j * (Constants.CARBON13 - Constants.CARBON));
        //            }
        //        }
        //    }
        //    else if (Form1.NEUCODE_SIXPLEX_LEU && Form1.CORRECTLEUDLOSS && countResidues('L', sequence) > 0)
        //    {
        //        // Fill first row
        //        uncorrectedMonoMass = (theoMasses[4, 0] + theoMasses[5, 0]) / 2.0;                

        //        // Fill rest of the array
        //        double mass;
        //        for (int i = 0; i < numLabels; i++)
        //        {
        //            mass = uncorrectedMonoMass - (i * (Constants.DEUTERIUM - Constants.HYDROGEN));
        //            for (int j = 0; j < numIsotopes; j++)
        //            {
        //                correctedMasses[i, j] = mass + (j * (Constants.CARBON13 - Constants.CARBON));
        //            }
        //        }
        //    }
        //    else if (Form1.NEUCODE_SIXPLEX_ARG && Form1.CORRECTARGPROLINE && countResidues('R', sequence) > 0 && countResidues('P', sequence) > 0)
        //    {
        //        double prolineConversionMass = 5 * (Constants.CARBON13 - Constants.CARBON);

        //        // Fill first row
        //        uncorrectedMonoMass = (theoMasses[4, 0] + theoMasses[5, 0]) / 2.0;

        //        // Fill rest of the array
        //        double mass;
        //        for (int i = 0; i < numLabels; i++)
        //        {
        //            mass = uncorrectedMonoMass + (i * prolineConversionMass);
        //            for (int j = 0; j < numIsotopes; j++)
        //            {
        //                correctedMasses[i, j] = mass + (j * (Constants.CARBON13 - Constants.CARBON));
        //            }
        //        }
        //    }

        //    // Find largest peaks
        //    ILabeledPeak peak;
        //    for (int i = 0; i < numLabels; i++)
        //    {
        //        for (int j = 0; j < numIsotopes; j++)
        //        {
        //            peak = largestPeak(correctedMasses[i, j], current, firstSearchMassRange, rawFile);
        //            if (peak != null)
        //            {
        //                correctedIntensities[i, j] = peak.Y;
        //            }
        //        }
        //    }

        //    // Find best peak combination
        //    double bestSummedIntensity = 0;
        //    int bestIndex = 0;
        //    double sum;
        //    for (int i = 0; i < numLabels; i++)
        //    {
        //        sum = 0;
        //        for (int j = 0; j < numIsotopes; j++)
        //        {
        //            sum += correctedIntensities[i, j];
        //        }
        //        if (sum > bestSummedIntensity)
        //        {
        //            bestSummedIntensity = sum;
        //            bestIndex = i;
        //        }
        //    }
        //    bestMonoMass = bestIndex;
        //    conversionFactor = bestMonoMass;

        //    return bestMonoMass;
        //}

        public bool GetTheoreticalResolvability(int cluster)
        {
            if (numLabels == 0) return false;

            else if (numIsotopologues < 2) return true;

            else
            {
                int charge = bestPSM.Charge;
                double experimentalSeparation = (spacingMassRange[1, 0].Mean) / (double)charge;
                int clusterIndex = cluster * numIsotopologues;
                double coefficient = (Math.Sqrt(2 * Math.Log(100.0 / QUANTSEPARATION))) / (Math.Sqrt(2 * Math.Log(2)));
                double theoreticalSeparation = coefficient * ((Mass.MzFromMass(theoMasses[clusterIndex, 0], charge) / (QUANTRESOLUTION * Math.Sqrt(400 / Mass.MzFromMass(theoMasses[clusterIndex, 0], charge)))));

                if (experimentalSeparation > theoreticalSeparation)
                {
                    return true;
                }
                return false;
            }
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
                            if (pair.complete[c, j])
                            {
                                pairTotalCount++;
                            }
                            // Pairs with fewer than half the total number of isotopologues
                            else if (Form1.NOISEBANDCAP && pair.peakCount[c, j] < peaksNeeded)
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
            if (numIsotopologues > 1 && pair.peakCount[cluster, isotope] != numIsotopologues)
            {
                return false;
            }
            else if (numIsotopologues < 2 && pair.peakCount[cluster, isotope] != numChannels)
            {
                return false;
            }

            // If there are more than 2 complete pairs, use their max intensity
            if (countCompleteIsotopes[cluster] > peaksNeeded)
            {
                max = maxCompleteIntensity;
            }
            else
            {
                max = maxIntensity;
            }

            if (numIsotopologues > 1)
            {
                channelIndex = cluster * numIsotopologues;
                int count;

                // For complete pairs
                if (pair.complete[cluster, isotope])
                {
                    count = 0;
                    for (int i = channelIndex; i < channelIndex + numIsotopologues; i++)
                    {
                        double intensityThreshold = INTENSITYCUTOFF * max[i, 0];
                        double peakIntensity = pair.peaks[i, isotope].GetDenormalizedIntensity(pair.injectionTime);
                        if (peakIntensity < intensityThreshold)
                        {
                            count++;
                        }
                    }
                }
                
                // For incomplete pairs, only consider non noise-band capped channels for quantitative filtering
                else
                {
                    count = 0;
                    for (int i = channelIndex; i < channelIndex + numIsotopologues; i++)
                    {
                        double intensityThreshold = INTENSITYCUTOFF * max[i, 0];
                        double peakIntensity = pair.peaks[i, isotope].GetDenormalizedIntensity(pair.injectionTime);
                        //if (pair.peakCount[cluster, isotope] > peaksNeeded && peakIntensity < intensityThreshold)
                        //{
                        //    count++;
                        //}
                        
                        if (pair.peaks[i, isotope].X > 0 && peakIntensity < intensityThreshold)
                        {
                            count++;
                        }
                    }
                }

                if (count > 0) return false;
                else return true;
            }
            else
            {
                int count;

                for (int i = 0; i < numChannels; i++)
                {
                    count = 0;
                    // For complete pairs
                    if (pair.complete[cluster, isotope])
                    {
                        double intensityThreshold = INTENSITYCUTOFF * max[i, 0];
                        double peakIntensity = pair.peaks[i, isotope].GetDenormalizedIntensity(pair.injectionTime);
                        if (peakIntensity < intensityThreshold)
                        {
                            count++;
                        }
                    }
                    // For incomplete pairs, only consider non noise-band capped channels for quantitative filtering
                    else
                    {
                        double intensityThreshold = INTENSITYCUTOFF * max[i, 0];
                        double peakIntensity = pair.peaks[i, isotope].GetDenormalizedIntensity(pair.injectionTime);
                        //if (pair.peakCount[cluster, isotope] > peaksNeeded && peakIntensity < intensityThreshold)
                        //{
                        //    count++;
                        //}

                        if (pair.peaks[i, isotope].X > 0 && peakIntensity < intensityThreshold)
                        {
                            count++;
                        }
                    }

                    if (count > 0) return false;
                    else return true;
                }
                
                //for (int i = 0; i < numChannels; i++)
                //{
                //    // For complete pairs
                //    if (pair.complete[cluster, isotope])
                //    {
                //        if (pair.peaks[i, isotope].GetDenormalizedIntensity(pair.injectionTime) < INTENSITYCUTOFF * max[i, 0])
                //        {
                //            return false;
                //        }
                //    }
                //    // For incomplete pairs, only consider non noise-band capped channels for quantitative filtering
                //    else
                //    {
                //        if (pair.peaks[i, isotope].X > 0 && (pair.peaks[i, isotope].GetDenormalizedIntensity(pair.injectionTime) < INTENSITYCUTOFF * max[i, 0]))
                //        {
                //            return false;
                //        }
                //    }
                //}
            }
            
            //// For complete pairs
            //if (complete)
            //{
            //    max = maxCompleteIntensity;

            //    if (Form1.NEUCODE)
            //    {
            //        channelIndex = cluster * numIsotopologues;
            //        for (int i = channelIndex; i < channelIndex + numIsotopologues; i++)
            //        {
            //            double intensityThreshold = INTENSITYCUTOFF * max[i, 0];
            //            double denormalizedIntensity = pair.peaks[i, isotope].GetDenormalizedIntensity(pair.injectionTime);
            //            if (denormalizedIntensity < intensityThreshold)
            //            {
            //                return false;
            //            }
            //        }
            //    }
            //    else
            //    {
            //        for (int i = 0; i < numChannels; i++)
            //        {
            //            double intensityThreshold = INTENSITYCUTOFF * max[i, 0];
            //            double denormalizedIntensity = pair.peaks[i, isotope].GetDenormalizedIntensity(pair.injectionTime);
            //            if (denormalizedIntensity < intensityThreshold)
            //            {
            //                return false;
            //            }
            //        }
            //    }   
            //}
            //// For incomplete pairs
            //else
            //{
            //    // Use maximum of complete pairs if more than 2 clusters are complete
            //    max = maxIntensity;

            //    if (Form1.NEUCODE)
            //    {
            //        channelIndex = cluster * numIsotopologues;
            //        for (int i = channelIndex; i < channelIndex + numIsotopologues; i++)
            //        {
            //            // For complete pairs
            //            if (pair.complete[cluster, isotope])
            //            {
            //                if (pair.peaks[i, isotope].GetDenormalizedIntensity(pair.injectionTime) < INTENSITYCUTOFF * max[i, 0])
            //                {
            //                    return false;
            //                }
            //            }
            //            // For incomplete pairs, only consider non noise-band capped channels for quantitative filtering
            //            if (!pair.complete[cluster, isotope])
            //            {
            //                if (pair.peaks[i, isotope].X > 0 && pair.peaks[i, isotope].GetDenormalizedIntensity(pair.injectionTime) < INTENSITYCUTOFF * max[i, 0])
            //                {
            //                    return false;
            //                }
            //            }
            //        }
            //    }
            //    else
            //    {
            //        for (int i = 0; i < numChannels; i++)
            //        {
            //            // For complete pairs
            //            if (pair.complete[cluster, isotope])
            //            {
            //                if (pair.peaks[i, isotope].GetDenormalizedIntensity(pair.injectionTime) < INTENSITYCUTOFF * max[i, 0])
            //                {
            //                    return false;
            //                }
            //            }
            //            // For incomplete pairs, only consider non noise-band capped channels for quantitative filtering
            //            if (!pair.complete[cluster, isotope])
            //            {
            //                if (pair.peaks[i, isotope].X > 0 && pair.peaks[i, isotope].GetDenormalizedIntensity(pair.injectionTime) < INTENSITYCUTOFF * max[i, 0])
            //                {
            //                    return false;
            //                }
            //            }
            //        }
            //    }
            //}
            return addIntensities;
        }

        /* Assembles all the peptide's quantitative information
         */
        public void quantify()
        {            
            if (numIsotopologues > 1)
            {
                int[,] final = new int[numClusters, 2];
                quantifiedNoiseIncluded = new bool[numClusters];
                finalQuantified = new int[numClusters];

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

                // First try to quantify based on complete sets of isotopologues
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

                // Purity correct all lysine-based NeuCode
                //if (Form1.LYSINEPURITYCORRECTION)
                //{
                //    Dictionary<int, double> purityCorrections = Form1.CORRECTIONFACTORS;
                //    // 2-plex NeuCode
                //    if (numIsotopologues == 2)
                //    {
                //        for (int i = 0; i < numChannels; i += numIsotopologues)
                //        {
                //            totalIntensity[i, numIsotopes] = totalIntensity[i, numIsotopes] / purityCorrections[1];
                //            totalIntensity[i + 1, numIsotopes] = totalIntensity[i + 1, numIsotopes] / purityCorrections[6];

                //            completeTotalIntensity[i, numIsotopes] = completeTotalIntensity[i, numIsotopes] / purityCorrections[1];
                //            completeTotalIntensity[i + 1, numIsotopes] = completeTotalIntensity[i + 1, numIsotopes] / purityCorrections[6];
                //        }
                //    }

                //    // 3-plex NeuCode
                //    else if (numIsotopologues == 3)
                //    {
                //        for (int i = 0; i < numChannels; i += numIsotopologues)
                //        {
                //            totalIntensity[i, numIsotopes] = totalIntensity[i, numIsotopes] / purityCorrections[1];
                //            totalIntensity[i + 1, numIsotopes] = totalIntensity[i + 1, numIsotopes] / purityCorrections[4];
                //            totalIntensity[i + 2, numIsotopes] = totalIntensity[i + 2, numIsotopes] / purityCorrections[6];

                //            completeTotalIntensity[i, numIsotopes] = completeTotalIntensity[i, numIsotopes] / purityCorrections[1];
                //            completeTotalIntensity[i + 1, numIsotopes] = completeTotalIntensity[i + 1, numIsotopes] / purityCorrections[4];
                //            completeTotalIntensity[i + 2, numIsotopes] = completeTotalIntensity[i + 2, numIsotopes] / purityCorrections[6];
                //        }
                //    }

                //    // 4-plex NeuCode
                //    else if (numIsotopologues == 4)
                //    {
                //        for (int i = 0; i < numChannels; i += numIsotopologues)
                //        {
                //            totalIntensity[i, numIsotopes] = totalIntensity[i, numIsotopes] / purityCorrections[1];
                //            totalIntensity[i + 1, numIsotopes] = totalIntensity[i + 1, numIsotopes] / purityCorrections[3];
                //            totalIntensity[i + 2, numIsotopes] = totalIntensity[i + 2, numIsotopes] / purityCorrections[5];
                //            totalIntensity[i + 3, numIsotopes] = totalIntensity[i + 3, numIsotopes] / purityCorrections[6];

                //            completeTotalIntensity[i, numIsotopes] = completeTotalIntensity[i, numIsotopes] / purityCorrections[1];
                //            completeTotalIntensity[i + 1, numIsotopes] = completeTotalIntensity[i + 1, numIsotopes] / purityCorrections[3];
                //            completeTotalIntensity[i + 2, numIsotopes] = completeTotalIntensity[i + 2, numIsotopes] / purityCorrections[5];
                //            completeTotalIntensity[i + 3, numIsotopes] = completeTotalIntensity[i + 3, numIsotopes] / purityCorrections[6];
                //        }
                //    }

                //    // 6-plex NeuCode
                //    else if (numIsotopologues == 6)
                //    {
                //        for (int i = 0; i < numChannels; i += numIsotopologues)
                //        {
                //            totalIntensity[i, numIsotopes] = totalIntensity[i, numIsotopes] / purityCorrections[1];
                //            totalIntensity[i + 1, numIsotopes] = totalIntensity[i + 1, numIsotopes] / purityCorrections[2];
                //            totalIntensity[i + 2, numIsotopes] = totalIntensity[i + 2, numIsotopes] / purityCorrections[3];
                //            totalIntensity[i + 3, numIsotopes] = totalIntensity[i + 3, numIsotopes] / purityCorrections[4];
                //            totalIntensity[i + 4, numIsotopes] = totalIntensity[i + 4, numIsotopes] / purityCorrections[5];
                //            totalIntensity[i + 5, numIsotopes] = totalIntensity[i + 5, numIsotopes] / purityCorrections[6];

                //            completeTotalIntensity[i, numIsotopes] = completeTotalIntensity[i, numIsotopes] / purityCorrections[1];
                //            completeTotalIntensity[i + 1, numIsotopes] = completeTotalIntensity[i + 1, numIsotopes] / purityCorrections[2];
                //            completeTotalIntensity[i + 2, numIsotopes] = completeTotalIntensity[i + 2, numIsotopes] / purityCorrections[3];
                //            completeTotalIntensity[i + 3, numIsotopes] = completeTotalIntensity[i + 3, numIsotopes] / purityCorrections[4];
                //            completeTotalIntensity[i + 4, numIsotopes] = completeTotalIntensity[i + 4, numIsotopes] / purityCorrections[5];
                //            completeTotalIntensity[i + 5, numIsotopes] = completeTotalIntensity[i + 5, numIsotopes] / purityCorrections[6];
                //        }
                //    }
                //}

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

                        if (final[c, 1] >= minimumPostQFPairs)
                        {
                            finalQuantified[c] = final[c, 1];
                            quantifiedNoiseIncluded[c] = false;
                            noQuantReason = NonQuantifiableType.Quantified;

                            foreach (MSDataFile rawFile in completePairs.Keys)
                            {
                                allPairs[rawFile] = completePairs[rawFile];
                            }

                            // Use only complete pairs for quantification
                            for (int i = channelIndex; i < channelIndex + numIsotopologues; i++)
                            {
                                totalIntensity[i, numIsotopes] = completeTotalIntensity[i, numIsotopes];
                            }

                            int index1 = channelIndex;
                            int index2 = index1 + 1;
                            // Use only complete pairs for quantification
                            for (int n = index2; n < channelIndex + numIsotopologues; n++)
                            {
                                lightInt = completeTotalIntensity[index1, numIsotopes];
                                heavyInt = completeTotalIntensity[n, numIsotopes];

                                if (lightInt > 0 && heavyInt > 0)
                                {
                                    heavyToLightRatioSum[n - (c + 1), 0] = heavyInt / lightInt;
                                }
                                index2++;
                            }
                        }
                        else if (final[c, 0] + final[c, 1] >= minimumPostQFPairs)
                        {
                            finalQuantified[c] = final[c, 0] + final[c, 1];
                            quantifiedNoiseIncluded[c] = true;
                            noQuantReason = NonQuantifiableType.Quantified;

                            // Use only complete pairs for quantification
                            for (int i = channelIndex; i < channelIndex + numIsotopologues; i++)
                            {
                                totalIntensity[i, numIsotopes] += completeTotalIntensity[i, numIsotopes];
                            }

                            int index1 = channelIndex;
                            int index2 = index1 + 1;
                            // Use only complete pairs for quantification
                            for (int n = index2; n < channelIndex + numIsotopologues; n++)
                            {
                                lightInt = totalIntensity[index1, numIsotopes];
                                heavyInt = totalIntensity[n, numIsotopes];

                                if (lightInt > 0 && heavyInt > 0)
                                {
                                    heavyToLightRatioSum[n - (c + 1), 0] = heavyInt / lightInt;
                                }
                                index2++;
                            }
                        }
                        else
                        {
                            if (countAllIsotopes[c] + countCompleteIsotopes[c] >= minimumPostQFPairs)
                            {
                                noQuantReason = NonQuantifiableType.WrongElutionProfiles;
                            }
                            else if (countAllIsotopes[c] + countCompleteIsotopes[c] > 0)
                            {
                                noQuantReason = NonQuantifiableType.NotEnoughMeasurements;
                            }
                            
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
            // Traditional SILAC
            else
            {
                int[,] final = new int[1, 2];
                quantifiedNoiseIncluded = new bool[1];
                finalQuantified = new int[1];

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

                // First try to quantify based on complete sets of isotopologues
                if (completePairs != null && completePairs.Count > 0)
                {
                    foreach (List<Pair> pairs in completePairs.Values)
                    {
                        foreach (Pair pair in pairs)
                        {
                            for (int j = 0; j < numIsotopes; j++)
                            {
                                if (Form1.QUANTFILTER)
                                {
                                    // Eliminate low-level pairs by quantitative filtering
                                    if (quantFilter(pair, j, 0, true))
                                    {
                                        for (int i = 0; i < numChannels; i++)
                                        {
                                            double denormalizedIntensity = pair.peaks[i, j].GetDenormalizedIntensity(pair.injectionTime);
                                            completeTotalIntensity[i, j] += denormalizedIntensity;
                                            completeTotalIntensity[i, numIsotopes] += denormalizedIntensity;
                                        }
                                        final[0, 1]++;
                                    }
                                }
                                else
                                {
                                    bool noNullPeaks = true;
                                    for (int i = 0; i < numChannels; i++)
                                    {
                                        if (pair.peaks[i, j] == null)
                                        {
                                            noNullPeaks = false;
                                        }
                                    }

                                    if (noNullPeaks)
                                    {
                                        for (int i = 0; i < numChannels; i++)
                                        {
                                            double denormalizedIntensity = pair.peaks[i, j].GetDenormalizedIntensity(pair.injectionTime);
                                            completeTotalIntensity[i, j] += denormalizedIntensity;
                                            completeTotalIntensity[i, numIsotopes] += denormalizedIntensity;
                                        }
                                        final[0, 1]++;
                                    }
                                }
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
                                if (Form1.QUANTFILTER)
                                {
                                    //Use a peak's intensity if it is not noise-band capped and its intensity is greater than 1/2e of the maximum intensity
                                    if (quantFilter(pair, j, 0, false))
                                    {
                                        for (int i = 0; i < numChannels; i++)
                                        {
                                            double denormalizedIntensity = pair.peaks[i, j].GetDenormalizedIntensity(pair.injectionTime);
                                            totalIntensity[i, j] += denormalizedIntensity;
                                            totalIntensity[i, numIsotopes] += denormalizedIntensity;

                                        }
                                        final[0, 0]++;
                                    }
                                }
                                else
                                {
                                    bool noNullPeaks = true;
                                    int realPeaks = 0;
                                    for (int i = 0; i < numChannels; i++)
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
                                        for (int i = 0; i < numChannels; i++)
                                        {
                                            if (pair.peaks[i, j] != null)
                                            {
                                                double denormalizedIntensity = pair.peaks[i, j].GetDenormalizedIntensity(pair.injectionTime);
                                                totalIntensity[i, j] += denormalizedIntensity;
                                                totalIntensity[i, numIsotopes] += denormalizedIntensity;
                                            }
                                        }
                                        final[0, 0]++;
                                    }
                                }
                            }
                        }
                    }
                }

                // Purity correct all lysine-based NeuCode
                //if (Form1.LYSINEPURITYCORRECTION)
                //{
                //    Dictionary<int, double> purityCorrections = Form1.CORRECTIONFACTORS;
                //    // 2-plex NeuCode
                //    if (numIsotopologues == 2)
                //    {
                //        for (int i = 0; i < numChannels; i += numIsotopologues)
                //        {
                //            totalIntensity[i, numIsotopes] = totalIntensity[i, numIsotopes] / purityCorrections[1];
                //            totalIntensity[i + 1, numIsotopes] = totalIntensity[i + 1, numIsotopes] / purityCorrections[6];

                //            completeTotalIntensity[i, numIsotopes] = completeTotalIntensity[i, numIsotopes] / purityCorrections[1];
                //            completeTotalIntensity[i + 1, numIsotopes] = completeTotalIntensity[i + 1, numIsotopes] / purityCorrections[6];
                //        }
                //    }

                //    // 3-plex NeuCode
                //    else if (numIsotopologues == 3)
                //    {
                //        for (int i = 0; i < numChannels; i += numIsotopologues)
                //        {
                //            totalIntensity[i, numIsotopes] = totalIntensity[i, numIsotopes] / purityCorrections[1];
                //            totalIntensity[i + 1, numIsotopes] = totalIntensity[i + 1, numIsotopes] / purityCorrections[4];
                //            totalIntensity[i + 2, numIsotopes] = totalIntensity[i + 2, numIsotopes] / purityCorrections[6];

                //            completeTotalIntensity[i, numIsotopes] = completeTotalIntensity[i, numIsotopes] / purityCorrections[1];
                //            completeTotalIntensity[i + 1, numIsotopes] = completeTotalIntensity[i + 1, numIsotopes] / purityCorrections[4];
                //            completeTotalIntensity[i + 2, numIsotopes] = completeTotalIntensity[i + 2, numIsotopes] / purityCorrections[6];
                //        }
                //    }

                //    // 4-plex NeuCode
                //    else if (numIsotopologues == 4)
                //    {
                //        for (int i = 0; i < numChannels; i += numIsotopologues)
                //        {
                //            totalIntensity[i, numIsotopes] = totalIntensity[i, numIsotopes] / purityCorrections[1];
                //            totalIntensity[i + 1, numIsotopes] = totalIntensity[i + 1, numIsotopes] / purityCorrections[3];
                //            totalIntensity[i + 2, numIsotopes] = totalIntensity[i + 2, numIsotopes] / purityCorrections[5];
                //            totalIntensity[i + 3, numIsotopes] = totalIntensity[i + 3, numIsotopes] / purityCorrections[6];

                //            completeTotalIntensity[i, numIsotopes] = completeTotalIntensity[i, numIsotopes] / purityCorrections[1];
                //            completeTotalIntensity[i + 1, numIsotopes] = completeTotalIntensity[i + 1, numIsotopes] / purityCorrections[3];
                //            completeTotalIntensity[i + 2, numIsotopes] = completeTotalIntensity[i + 2, numIsotopes] / purityCorrections[5];
                //            completeTotalIntensity[i + 3, numIsotopes] = completeTotalIntensity[i + 3, numIsotopes] / purityCorrections[6];
                //        }
                //    }

                //    // 6-plex NeuCode
                //    else if (numIsotopologues == 6)
                //    {
                //        for (int i = 0; i < numChannels; i += numIsotopologues)
                //        {
                //            totalIntensity[i, numIsotopes] = totalIntensity[i, numIsotopes] / purityCorrections[1];
                //            totalIntensity[i + 1, numIsotopes] = totalIntensity[i + 1, numIsotopes] / purityCorrections[2];
                //            totalIntensity[i + 2, numIsotopes] = totalIntensity[i + 2, numIsotopes] / purityCorrections[3];
                //            totalIntensity[i + 3, numIsotopes] = totalIntensity[i + 3, numIsotopes] / purityCorrections[4];
                //            totalIntensity[i + 4, numIsotopes] = totalIntensity[i + 4, numIsotopes] / purityCorrections[5];
                //            totalIntensity[i + 5, numIsotopes] = totalIntensity[i + 5, numIsotopes] / purityCorrections[6];

                //            completeTotalIntensity[i, numIsotopes] = completeTotalIntensity[i, numIsotopes] / purityCorrections[1];
                //            completeTotalIntensity[i + 1, numIsotopes] = completeTotalIntensity[i + 1, numIsotopes] / purityCorrections[2];
                //            completeTotalIntensity[i + 2, numIsotopes] = completeTotalIntensity[i + 2, numIsotopes] / purityCorrections[3];
                //            completeTotalIntensity[i + 3, numIsotopes] = completeTotalIntensity[i + 3, numIsotopes] / purityCorrections[4];
                //            completeTotalIntensity[i + 4, numIsotopes] = completeTotalIntensity[i + 4, numIsotopes] / purityCorrections[5];
                //            completeTotalIntensity[i + 5, numIsotopes] = completeTotalIntensity[i + 5, numIsotopes] / purityCorrections[6];
                //        }
                //    }
                //}

                if (final[0, 1] >= minimumPostQFPairs)
                {
                    finalQuantified[0] = final[0, 1];
                    quantifiedNoiseIncluded[0] = false;
                    noQuantReason = NonQuantifiableType.Quantified;

                    totalIntensity = completeTotalIntensity;
                    
                    int index1 = 0;
                    int index2 = index1 + 1;
                    
                    // Use only complete pairs for quantification
                    for (int i = index2; i < numChannels - 1; i++)
                    {
                        // Use only complete pairs for quantification
                        lightInt = completeTotalIntensity[index1, numIsotopes];
                        heavyInt = completeTotalIntensity[i, numIsotopes];

                        if (lightInt > 0 && heavyInt > 0)
                        {
                            heavyToLightRatioSum[index2 - 1, 0] = heavyInt / lightInt;
                        }
                    }
                }
                else if (final[0, 0] + final[0, 1] >= minimumPostQFPairs)
                {
                    finalQuantified[0] = final[0, 0] + final[0, 1];
                    quantifiedNoiseIncluded[0] = true;
                    noQuantReason = NonQuantifiableType.Quantified;
                    
                    int index1 = 0;
                    int index2 = index1 + 1;
                    totalIntensity[index1, numIsotopes] += completeTotalIntensity[index1, numIsotopes];
                    
                    // Use only complete pairs for quantification
                    for (int i = index2; i < numChannels - 1; i++)
                    {
                        totalIntensity[i, numIsotopes] += completeTotalIntensity[i, numIsotopes];
                        // Use only complete pairs for quantification
                        lightInt = totalIntensity[index1, numIsotopes];
                        heavyInt = totalIntensity[i, numIsotopes];

                        if (lightInt > 0 && heavyInt > 0)
                        {
                            heavyToLightRatioSum[index2 - 1, 0] = heavyInt / lightInt;
                        }
                    }
                }
                else
                {
                    if (countAllIsotopes[0] + countCompleteIsotopes[0] >= minimumPostQFPairs)
                    {
                        noQuantReason = NonQuantifiableType.WrongElutionProfiles;
                    }
                    else if (countAllIsotopes[0] + countCompleteIsotopes[0] > 0)
                    {
                        noQuantReason = NonQuantifiableType.NotEnoughMeasurements;
                    }

                    // Not able to quantify
                    for (int i = 0; i < numChannels; i++)
                    {
                        completeTotalIntensity[i, numIsotopes] = 0;
                        totalIntensity[i, numIsotopes] = 0;
                    }
                }

                //int[,] final = new int[1, 2];
                //quantifiedNoiseIncluded = new bool[1];
                //finalQuantified = new int[1];

                //if (completePairs != null && completePairs.Count > 0)
                //{
                //    foreach (List<Pair> pairs in completePairs.Values)
                //    {
                //        foreach (Pair pair in pairs)
                //        {
                //            for (int j = 0; j < numIsotopes; j++)
                //            {
                //                if (Form1.QUANTFILTER)
                //                {
                //                    // Eliminate low-level pairs by quantitative filtering
                //                    if (quantFilter(pair, j, 0, true))
                //                    {
                //                        for (int i = 0; i < numChannels; i++)
                //                        {
                //                            double denormalizedIntensity = pair.peaks[i, j].GetDenormalizedIntensity(pair.injectionTime);
                //                            completeTotalIntensity[i, j] += denormalizedIntensity;
                //                            completeTotalIntensity[i, numIsotopes] += denormalizedIntensity;
                //                        }
                //                        final[0, 1]++;
                //                    }
                //                }
                //                else
                //                {
                //                    bool noNullPeaks = true;
                //                    for (int i = 0; i < numChannels; i++)
                //                    {
                //                        if (pair.peaks[i, j] == null)
                //                        {
                //                            noNullPeaks = false;
                //                        }
                //                    }

                //                    if (noNullPeaks)
                //                    {
                //                        for (int i = 0; i < numChannels; i++)
                //                        {
                //                            double denormalizedIntensity = pair.peaks[i, j].GetDenormalizedIntensity(pair.injectionTime);
                //                            completeTotalIntensity[i, j] += denormalizedIntensity;
                //                            completeTotalIntensity[i, numIsotopes] += denormalizedIntensity;
                //                        }
                //                        final[0, 1]++;
                //                    }
                //                }
                //            } // End isotope loop
                //        } // End pair loop
                //    } // End pair list loop
                //}
                //if (allPairs != null && allPairs.Count > 0)
                //{
                //    foreach (List<Pair> pairs in allPairs.Values)
                //    {
                //        foreach (Pair pair in pairs)
                //        {
                //            for (int j = 0; j < numIsotopes; j++)
                //            {
                //                if (Form1.QUANTFILTER)
                //                {
                //                    //Use a peak's intensity if it is not noise-band capped and its intensity is greater than 1/2e of the maximum intensity

                //                    if (quantFilter(pair, j, 0, false))
                //                    {
                //                        for (int i = 0; i < numChannels; i++)
                //                        {
                //                            double denormalizedIntensity = pair.peaks[i, j].GetDenormalizedIntensity(pair.injectionTime);
                //                            totalIntensity[i, j] += denormalizedIntensity;
                //                            totalIntensity[i, numIsotopes] += denormalizedIntensity;

                //                        }
                //                        final[0, 0]++;
                //                    }
                //                }
                //                else
                //                {
                //                    bool noNullPeaks = true;
                //                    int realPeaks = 0;
                //                    for (int i = 0; i < numChannels; i++)
                //                    {
                //                        if (pair.peaks[i, j] == null)
                //                        {
                //                            noNullPeaks = false;
                //                        }
                //                        else
                //                        {
                //                            realPeaks++;
                //                        }
                //                    }

                //                    if (noNullPeaks)
                //                    {
                //                        for (int i = 0; i < numChannels; i++)
                //                        {
                //                            if (pair.peaks[i, j] != null)
                //                            {
                //                                double denormalizedIntensity = pair.peaks[i, j].GetDenormalizedIntensity(pair.injectionTime);
                //                                totalIntensity[i, j] += denormalizedIntensity;
                //                                totalIntensity[i, numIsotopes] += denormalizedIntensity;
                //                            }
                //                        }
                //                        final[0, 0]++;
                //                    }
                //                }
                //            }
                //        }
                //    }
                //}

                //double lightInt;
                //double heavyInt;
                //bool quantified = false;
                //int minimumTotalPairs = 3;
                //int minimumNoNBCPairs = 3;
                //int minimumPostQFPairs = 3;

                //if (final[0, 1] >= minimumPostQFPairs)
                //{
                //    finalQuantified[0] = final[0, 1];
                //    quantifiedNoiseIncluded[0] = false;

                //    for (int i = 1; i < numChannels; i++)
                //    {
                //        if (Form1.NOISEBANDCAP)
                //        {
                //            totalIntensity = completeTotalIntensity;
                //        }
                //        // Use only complete pairs for quantification
                //        int index1 = 0;
                //        int index2 = i;

                //        lightInt = completeTotalIntensity[index1, numIsotopes];
                //        heavyInt = completeTotalIntensity[index2, numIsotopes];

                //        if (lightInt > 0 && heavyInt > 0)
                //        {
                //            heavyToLightRatioSum[i - 1, 0] = heavyInt / lightInt;
                //        }
                //    }
                //}
                //else if (final[0, 0] >= minimumPostQFPairs)
                //{
                //    finalQuantified[0] = final[0, 0];
                //    quantifiedNoiseIncluded[0] = true;

                //    for (int i = 1; i < numChannels; i++)
                //    {
                //        int index1 = 0;
                //        int index2 = i;

                //        lightInt = totalIntensity[index1, numIsotopes];
                //        heavyInt = totalIntensity[index2, numIsotopes];

                //        if (lightInt > 0 && heavyInt > 0)
                //        {
                //            heavyToLightRatioSum[i - 1, 0] = heavyInt / lightInt;
                //        }
                //    }
                //}
                //else
                //{
                //    // Not able to quantify
                //    for (int i = 0; i < numChannels; i++)
                //    {
                //        completeTotalIntensity[i, numIsotopes] = 0;
                //        totalIntensity[i, numIsotopes] = 0;
                //    }
                //}
            }  
        }

        public static double calculateAverage(List<double> numbers)
        {
            double sum = 0;
            double count = (double)numbers.Count;

            if (count == 0) return 0;

            foreach (double number in numbers)
            {
                sum += number;
            }

            return sum / count;
        }

        public static double calculateMedian(List<double> numbers)
        {
            int count = numbers.Count;

            numbers.Sort();

            if (count == 0)
            {
                return 0;
            }
            else if (count % 2 != 0)
            {
                return numbers[count / 2];
            }
            else
            {
                return (numbers[count / 2] + numbers[(count / 2) - 1]) / 2.0;
            }
        }
    }
}
