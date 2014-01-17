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
    public class PeptideID
    {
        // Constants
        //public static MassTolerance SILAC = new MassTolerance(MassToleranceType.PPM, 20.0); // Individual SILAC peak tolerance (for lower resolution MS1 scans)
        //public static MassTolerance NEUCODE = new MassTolerance(MassToleranceType.PPM, 10.0); // Individual NeuCode peak tolerance (for higher resolution MS1 scans)
        public static double INTENSITYCUTOFF = 1.0 / Math.E; // Intensity threshold for quantitation filtering to eliminate low-level peaks
        public static double SPACINGPERCENTAGEERROR = 0.20;

        // User-specified tolerances
        public double minSignalToNoise
        {
            get
            {
                return Form1.MINIMUMSN;
            }
        }
        public double maxRawIntensity
        {
            get
            {
                if (Form1.MAXIMUMNL > 0)
                {
                    return Form1.MAXIMUMNL * 1000000.0;
                }
                else
                {
                    return double.PositiveInfinity;
                }                
            }
        }
        public double quantResolution
        {
            get
            {
                if (Form1.FUSION)
                {
                    return Form1.QUANTRESOLUTION * 1000.0 / Math.Sqrt(2.0);
                }
                else
                {
                    return Form1.QUANTRESOLUTION * 1000.0;
                }
            }
        }            
        public double quantSeparation 
        {
            get
            {
                return Form1.THEORETICALSEPARATION;
            }
        }
            
        public MassTolerance ppmTolerance 
        {
            get
            {
                return new MassTolerance(MassToleranceType.PPM, Form1.TOLERANCE);
            }
        }
            
        // Class members
        // Experimental design information
        public int numChannels; // # of total channels in the experiment (# clusters * # isotopologues)
        public int numIsotopes; // # of isotopes to include
        public int numClusters; // # of clusters
        public int numIsotopologues; // # of isotopologues within each cluster
        public int numPeptideVersions { get; set; }
        
        // Identification information
        public Dictionary<MSDataFile, List<PeptideSpectralMatch>> PSMs; // Organizes a peptide's PSMs based on raw file
        public PeptideSpectralMatch PSM; // A peptide's first PSM
        public Dictionary<MSDataFile, PeptideSpectralMatch> bestPSMs;
        public PeptideSpectralMatch bestPSM;
        //{
        //    get
        //    {
        //        PeptideSpectralMatch best = null;
        //        if (bestPSMs.Count > 0)
        //        {
        //            foreach (PeptideSpectralMatch psm in bestPSMs.Values)
        //            {
        //                if (best == null || psm.EValue < best.EValue) best = psm;
        //            }

        //            //if (numIsotopologues > 1) best = maximizeResolvability();
        //        }
        //        return best;
        //    }
        //}
        //public int numIsotopologueLabels; // The number of quantitative labels carried by the peptide
        //public int numClusterLabels;
        //public int numLabels
        //{
        //    get
        //    {
        //        if (numIsotopologues > 1) return numIsotopologueLabels;
        //        else return numClusterLabels;
        //    }
        //}
        public bool clusterLabel;
        public bool isotopologueLabel;
        public bool labeled;
        public string sequence; // A peptide's sequence
        public string sequenceNoMods;
        public string sequenceItoL;
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

                scanNumbers = "\"" + scanNumbers + "\"";
                return scanNumbers.ToString();
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
                            rawFiles += PSMs.ElementAt(i).Key + ":";
                        }
                        else
                        {
                            rawFiles += PSMs.ElementAt(i).Key;
                        }
                    }
                }
                rawFiles = "\"" + rawFiles + "\"";
                return rawFiles.ToString();
            }
        } // All the raw files that produced PSMs
        public string chargeStates
        {
            get
            {
                string chargeStates = "";
                if (PSMs != null && PSMs.Count > 0)
                {
                    foreach (MSDataFile rawFile in PSMs.Keys) //For every raw file
                    {
                        List<PeptideSpectralMatch> psms = PSMs[rawFile];
                        foreach (PeptideSpectralMatch PSM in psms) //For every PSM
                        {
                            if (psms.ElementAt(psms.Count - 1) == PSM)
                            {
                                chargeStates += PSM.Charge;
                            }
                            else
                            {
                                chargeStates += PSM.Charge + ",";
                            }
                        } //End for every PSM

                        if (PSMs.Keys.ElementAt(PSMs.Keys.Count - 1) != rawFile)
                        {
                            chargeStates += ";";
                        }
                    } //End for every raw file
                }
                chargeStates = "\"" + chargeStates + "\"";
                return chargeStates.ToString();
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
                    for (int i = 0; i < numPeptideVersions; i++)
                    {
                        for (int j = 0; j < numIsotopes; j++)
                        {
                            if (theoMasses[i, 0] > 0)
                            {
                                adjustedMass = ((theoMasses[i, 0] * ppmError) / 1000000.0) + theoMasses[i, 0];
                                adjustedPrecursorMasses[i, j] = adjustedMass + (double)j * (Constants.Carbon13 - Constants.Carbon);
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
                    for (int i = 0; i < numPeptideVersions; i++)
                    {
                        for (int j = 0; j < numIsotopes; j++)
                        {
                            adjustedMass = theoMasses[i, 0];
                            adjustedPrecursorMasses[i, j] = adjustedMass + (double)j * (Constants.Carbon13 - Constants.Carbon);
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
                    if (numIsotopologues > 1 && isotopologueLabel)
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
                            if (numChannels > numPeptideVersions && i >= numPeptideVersions)
                            {
                                averageTheoMZ[i] = 0.0;
                            }
                            else
                            {
                                averageTheoMZ[i] = Mass.MzFromMass(theoMasses[i, 0], charge);
                            }
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
                            if (averageTheoMZ[c] > 0)
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
                        }
                        return averageAdjustedTheoMZ;
                    }
                    else
                    {
                        averageAdjustedTheoMZ = new double[numChannels];
                        for (int i = 0; i < numChannels; i++)
                        {
                            if (averageTheoMZ[i] > 0)
                            {
                                averageAdjustedTheoMZ[i] = Mass.MzFromMass(adjustedTheoMasses[i, 0], charge);
                            }
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
                        //if (numChannels == 3) numPeaks = 2;
                    }
                    else if (numIsotopologues > 1 && Form1.CROSSCLUSTERQUANT)
                    {
                        numPeaks = numChannels / 2;
                    }
                    else
                    {
                        numPeaks = numIsotopologues / 2;
                        //if (numIsotopologues == 3) numPeaks = 2;
                    }
                }
                else
                {
                    if (numIsotopologues < 2)
                    {
                        numPeaks = numChannels;
                    }
                    else if (numIsotopologues > 1 && Form1.CROSSCLUSTERQUANT)
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
                double[,] max = new double[numChannels, 2];
                List<Pair> pairsCombined = new List<Pair>();
                List<Pair> intensitySortedPairs = new List<Pair>();

                if (allPairs != null && allPairs.Count > 0)
                {
                    foreach (List<Pair> pairList in allPairs.Values)
                    {
                        foreach (Pair pair in pairList)
                        {
                            pairsCombined.Add(pair);
                        }
                    }
                }
                if (completePairs != null && completePairs.Count > 0)
                {
                    foreach (List<Pair> pairList in completePairs.Values)
                    {
                        foreach (Pair pair in pairList)
                        {
                            pairsCombined.Add(pair);
                        }
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
                            max[i, 1] = intensitySortedPairs[0].RetentionTime;
                            break;
                        case 2:
                            max[i, 0] = (intensitySortedPairs[0].GetMaxChannelIntensity(i) + intensitySortedPairs[1].GetMaxChannelIntensity(i)) / 2.0;
                            max[i, 1] = intensitySortedPairs[0].RetentionTime;
                            break;
                        case 3:
                            max[i, 0] = intensitySortedPairs[1].GetMaxChannelIntensity(i);
                            max[i, 1] = intensitySortedPairs[0].RetentionTime;
                            break;
                        case 4:
                            max[i, 0] = (intensitySortedPairs[1].GetMaxChannelIntensity(i) + intensitySortedPairs[2].GetMaxChannelIntensity(i)) / 2.0;
                            max[i, 1] = intensitySortedPairs[0].RetentionTime;
                            break;
                        case 5:
                            max[i, 0] = intensitySortedPairs[2].GetMaxChannelIntensity(i);
                            max[i, 1] = intensitySortedPairs[0].RetentionTime;
                            break;
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
                                max[i, 1] = intensitySortedPairs[0].RetentionTime;
                                break;
                            case 2:
                                max[i, 0] = (intensitySortedPairs[0].GetMaxChannelIntensity(i) + intensitySortedPairs[1].GetMaxChannelIntensity(i)) / 2.0;
                                max[i, 1] = intensitySortedPairs[0].RetentionTime;
                                break;
                            case 3:
                                max[i, 0] = intensitySortedPairs[1].GetMaxChannelIntensity(i);
                                max[i, 1] = intensitySortedPairs[0].RetentionTime;
                                break;
                            case 4:
                                max[i, 0] = (intensitySortedPairs[1].GetMaxChannelIntensity(i) + intensitySortedPairs[2].GetMaxChannelIntensity(i)) / 2.0;
                                max[i, 1] = intensitySortedPairs[0].RetentionTime;
                                break;
                            case 5:
                                max[i, 0] = intensitySortedPairs[2].GetMaxChannelIntensity(i);
                                max[i, 1] = intensitySortedPairs[0].RetentionTime;
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
                double[,] max = new double[numChannels, 2];
                List<Pair> pairsCombined = new List<Pair>();
                List<Pair> intensitySortedPairs = new List<Pair>();

                if (completePairs != null && completePairs.Count > 0)
                {
                    foreach (List<Pair> pairList in completePairs.Values)
                    {
                        foreach (Pair pair in pairList)
                        {
                            pairsCombined.Add(pair);
                        }
                    }
                }

                if (pairsCombined.Count == 0)
                {
                    if (allPairs != null && allPairs.Count > 0)
                    {
                        foreach (List<Pair> pairList in allPairs.Values)
                        {
                            foreach (Pair pair in pairList)
                            {
                                pairsCombined.Add(pair);
                            }
                        }
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
                            max[i, 1] = intensitySortedPairs[0].RetentionTime;
                            break;
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
                                max[i, 1] = intensitySortedPairs[0].RetentionTime;
                                break;
                        }
                    }
                }
                return max;
            }
        }
        public double[,] completeTotalIntensity; // Keeps track of each channel's total intensity considering only complete pairs
        public double[,] totalIntensity; // Keeps track of each channel's total intensity considering all pairs
        public double[] medianCompleteTotalIntensity;
        public double[] medianTotalIntensity;
        public MassTolerance firstSearchMassRange; // The mass range used to search for peaks for systematic error determination
        public double[] heavyToLightRatioSum; // Computes the quantitative ratio between channels
        public double[,] heavyToLightRatioMedian;
        public Dictionary<MSDataFile, List<Pair>> allPairs;
        public List<double>[,] channelMZValues
        {
            get
            {
                List<double>[,] MZ = new List<double>[numChannels, numIsotopes];

                for (int c = 0; c < numChannels; c++)
                {
                    for (int i = 0; i < numIsotopes; i++)
                    {
                        MZ[c, i] = new List<double>();
                    }
                }

                foreach (List<Pair> pairList in completePairs.Values)
                {
                    foreach (Pair pair in pairList)
                    {
                        foreach (IsotopePair isotopePair in pair.IsotopePairs)
                        {
                            if (isotopePair != null)
                            {
                                for (int c = 0; c < numChannels; c++)
                                {
                                    if (isotopePair.ChannelPeaks[c] != null)
                                    {
                                        MZ[c, isotopePair.Isotope].Add(isotopePair.ChannelPeaks[c].X);
                                    }
                                }
                            }
                        }
                    }
                }

                bool[,] completePairFound = new bool[numChannels, numIsotopes];

                for (int c = 0; c < numChannels; c++)
                {
                    for (int i = 0; i < numIsotopes; i++)
                    {
                        completePairFound[c, i] = false;
                        if (MZ[c, i].Count > 0) completePairFound[c, i] = true;
                    }
                }

                foreach (List<Pair> pairList in allPairs.Values)
                {
                    foreach (Pair pair in pairList)
                    {
                        foreach (IsotopePair isotopePair in pair.IsotopePairs)
                        {
                            if (isotopePair != null)
                            {
                                for (int c = 0; c < numChannels; c++)
                                {
                                    if (isotopePair.ChannelPeaks[c] != null && !completePairFound[c, isotopePair.Isotope])
                                    {
                                        MZ[c, isotopePair.Isotope].Add(isotopePair.ChannelPeaks[c].X);
                                    }
                                }
                            }
                        }
                    }
                }

                for (int c = 0; c < numChannels; c++)
                {
                    for (int i = 0; i < numIsotopes; i++)
                    {
                        MZ[c, i].Sort();
                    }
                }

                return MZ;
            }
        }
        public Chromatogram[,] channelXICs
        {
            get
            {
                Chromatogram[,] XICs = new Chromatogram[numChannels, numIsotopes];
                for (int c = 0; c < numChannels; c++)
                {
                    for (int i = 0; i < numIsotopes; i++)
                    {
                        XICs[c, i] = fullScanList.GetChromatogram(ChromatogramType.MzRange, channelXICRange[c, i]); 
                    }
                }
                return XICs;
            }
        }
        public double[,] averageChannelMZ
        {
            get
            {
                double[,] averageMZ = new double[numChannels, numIsotopes];
                for (int c = 0; c < numChannels; c++)
                {
                    for (int i = 0; i < numIsotopes; i++)
                    {
                        List<double> mZList = channelMZValues[c, i];
                        double total = 0.0;
                        int count = 0;
                        foreach (double MZ in mZList)
                        {
                            total += MZ;
                            count++;
                        }
                        averageMZ[c, i] = total / count;
                    }
                }
                return averageMZ;
            }
        }
        public MassTolerance[,] channelDeviationPPM
        {
            get
            {
                MassTolerance[,] stdDevPPM = new MassTolerance[numChannels, numIsotopes];
                for (int c = 0; c < numChannels; c++)
                {
                    for (int i = 0; i < numIsotopes; i++)
                    {
                        List<double> mZList = channelMZValues[c, i];
                        double averageMZ = averageChannelMZ[c, i];
                        double numerator = 0.0;
                        int denominator = mZList.Count - 1;
                        foreach (double MZ in mZList)
                        {
                            numerator += Math.Pow(MZ - averageMZ, 2);
                        }
                        double stdDev = Math.Sqrt(numerator / denominator);
                        stdDevPPM[c, i] = new MassTolerance(MassToleranceType.PPM, (stdDev * 2000000) / averageMZ);
                    }
                }
                return stdDevPPM;
            }
        }
        public MassRange[,] channelXICRange
        {
            get
            {
                MassRange[,] channelXICs = new MassRange[numChannels, numIsotopes];
                for (int c = 0; c < numChannels; c++)
                {
                    for (int i = 0; i < numIsotopes; i++)
                    {
                        if (averageChannelMZ[c, i] > 0)
                        {
                            if (channelDeviationPPM[c, i].Value > 0.0)
                            {
                                channelXICs[c, i] = new MassRange(averageChannelMZ[c, i], channelDeviationPPM[c, i]);
                            }
                            else
                            {
                                channelXICs[c, i] = new MassRange(averageChannelMZ[c, i], new MassTolerance(MassToleranceType.PPM, 2.0));
                            }
                        }
                        else
                        {
                            double adjustedTheoMZ = Mass.MzFromMass(adjustedTheoMasses[c, i], bestPSM.Charge);
                            channelXICs[c, i] = new MassRange(adjustedTheoMZ, new MassTolerance(MassToleranceType.PPM, 2.0));
                        }
                    }
                }
                return channelXICs;
            }
        }
        public int countAllPairs
        {
            get
            {
                int count = 0;
                if (allPairs != null && allPairs.Count > 0)
                {
                    foreach (MSDataFile rawFile in allPairs.Keys)
                    {
                        foreach (Pair pair in allPairs[rawFile])
                        {
                            if (pair != null && pair.IsotopePairs.Count > 0)
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
            get
            {
                int[] count = null;
                if (numIsotopologues > 1 && numClusters > 1)
                {
                    if (allPairs != null && allPairs.Count > 0)
                    {
                        count = new int[numClusters];
                        foreach (List<Pair> pairs in allPairs.Values)
                        {
                            foreach (Pair pair in pairs)
                            {
                                for (int c = 0; c < numClusters; c++)
                                {
                                    foreach (IsotopePair isotopePair in pair.IsotopePairs)
                                    {
                                        if (isotopePair.PeakCountByCluster[c] >= peaksNeeded)
                                        {
                                            count[c]++;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                else
                {
                    if (allPairs != null && allPairs.Count > 0)
                    {
                        count = new int[1];
                        foreach (List<Pair> pairs in allPairs.Values)
                        {
                            foreach (Pair pair in pairs)
                            {
                                foreach (IsotopePair isotopePair in pair.IsotopePairs)
                                {
                                    if (isotopePair.PeakCountByCluster[0] >= peaksNeeded)
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
            get
            {
                int count = 0;
                if (completePairs != null && completePairs.Count > 0)
                {
                    foreach (MSDataFile rawFile in completePairs.Keys)
                    {
                        foreach (Pair pair in completePairs[rawFile])
                        {
                            if (pair != null && pair.IsotopePairs.Count > 0)
                            {
                                count++;
                            }
                        }
                    }
                }
                return count;
            }
        }
        public int[] countCompleteIsotopes
        {
            get
            {
                int[] count = null;
                if (numIsotopologues > 1 && numClusters > 1)
                {
                    if (completePairs != null && completePairs.Count > 0)
                    {
                        count = new int[numClusters];
                        foreach (List<Pair> pairs in completePairs.Values)
                        {
                            foreach (Pair pair in pairs)
                            {
                                for (int c = 0; c < numClusters; c++)
                                {
                                    foreach (IsotopePair isotopePair in pair.IsotopePairs)
                                    {
                                        if (isotopePair.PeakCountByCluster[c] >= peaksNeeded)
                                        {
                                            count[c]++;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                else
                {
                    if (completePairs != null && completePairs.Count > 0)
                    {
                        count = new int[1];
                        foreach (List<Pair> pairs in completePairs.Values)
                        {
                            foreach (Pair pair in pairs)
                            {
                                foreach (IsotopePair isotopePair in pair.IsotopePairs)
                                {
                                    if (isotopePair.PeakCountByCluster[0] >= peaksNeeded)
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
        public bool[] quantifiedNoiseIncluded;
        public int[,] missingChannelPattern; 
        public int[] finalQuantified;
        public MassRange[,] spacingMassRange
        {
            get
            {
                MassRange[,] spacing;
                double theoSpacing;

                if (numIsotopologues > 1 && isotopologueLabel)
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
                                theoSpacing = theoMasses[r,0] - theoMasses[c,0];
                                if (theoSpacing < 0) theoSpacing = theoSpacing * -1.0;
                                spacing[r, c] = new MassRange(theoSpacing - (2 * (SPACINGPERCENTAGEERROR / 2) * theoSpacing), theoSpacing + (1 * (SPACINGPERCENTAGEERROR / 2) * theoSpacing));
                            }
                        }     
                    }
                    return spacing;
                }
                else if (numIsotopologues < 2 && clusterLabel)
                {
                    spacing = new MassRange[numChannels, numChannels];
                    theoSpacing = 0;

                    for (int r = 0; r < numChannels; r++)
                    {
                        for (int c = 0; c < numChannels; c++)
                        {
                            if (r == c) spacing[r, c] = null;
                            else
                            {
                                theoSpacing = theoMasses[r,0] - theoMasses[c,0];
                                if (theoSpacing < 0) theoSpacing = theoSpacing * -1.0;
                                spacing[r, c] = new MassRange(theoSpacing, new MassTolerance(MassToleranceType.DA, 0.01 * theoSpacing));
                            }
                            
                            
                            //if (r == c) spacing[r, c] = null;
                            //else if (r == 0)
                            //{
                            //    theoSpacing = 0 - ((double)numClusterLabels * Form1.CLUSTERLABELS[c - 1].MonoisotopicMass);
                            //    if (theoSpacing < 0) theoSpacing = theoSpacing * -1.0;
                            //    spacing[r, c] = new MassRange(theoSpacing, new MassTolerance(MassToleranceType.DA, 0.020 * numClusterLabels));
                            //}
                            //else if (c == 0)
                            //{
                            //    theoSpacing = ((double)numClusterLabels * Form1.CLUSTERLABELS[r - 1].MonoisotopicMass) - 0;
                            //    if (theoSpacing < 0) theoSpacing = theoSpacing * -1.0;
                            //    spacing[r, c] = new MassRange(theoSpacing, new MassTolerance(MassToleranceType.DA, 0.020 * numClusterLabels));
                            //}
                            //else
                            //{
                            //    theoSpacing = (double)numClusterLabels * (Form1.CLUSTERLABELS[r - 1].MonoisotopicMass - Form1.CLUSTERLABELS[c - 1].MonoisotopicMass);
                            //    if (theoSpacing < 0) theoSpacing = theoSpacing * -1.0;
                            //    spacing[r, c] = new MassRange(theoSpacing, new MassTolerance(MassToleranceType.DA, 0.020 * numClusterLabels));
                            //}
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
                    if (clusterLabel) return true;
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
        public bool allClustersResolvable
        {
            get
            {
                if (numIsotopologues < 2)
                {
                    if (clusterLabel) return true;
                    else return false;
                }
                else
                {
                    for (int i = 0; i < numClusters; i++)
                    {
                        if (!GetTheoreticalResolvability(i))
                        {
                            return false;
                        }
                    }
                    return true;
                }
            }
        }
        public int channelsDetected
        {
            get
            {
                if (noQuantReason == NonQuantifiableType.Quantified)
                {
                    if (quantifiedNoiseIncluded[0]) return 1;
                    else
                    {
                        return numChannels;
                    }
                }
                else return 0;
            }
        }
        public double[] theoreticalIsotopicDistribution;

        // Coalescence information
        public bool coalescenceDetected;
        public List<double> coalescedPeakIntensities;
        public Dictionary<int, List<double>> missingChannelsSN;

        // Creates a PeptideID based on a database search peptide identification
        public PeptideID(int scanNumber, int charge, double eValue, string sequenceOriginal, MSDataFile rawFile, string mods, string fileNameID)
        {
            string[] modifications = mods.Split(',');
            bool[] lysinePTM = new bool[modifications.Length];
            int[] lysineSites = new int[modifications.Length];
            if (Form1.ACETYLATION || Form1.UBIQUITIN && mods.Length > 1)
            {
                for (int m = 0; m < modifications.Length; m++)
                {
                    lysineSites[m] = int.Parse(modifications[m].Split(':')[1]);
                    if (modifications[m].ToLower().Contains("acetyl") || modifications[m].ToLower().Contains("ubiquit"))
                    {
                        lysinePTM[m] = true;
                    }
                }
            }
            
            // Local variables
            List<PeptideSpectralMatch> psms;
            Peptide peptide;
            Peptide[] peptideVersions = null;

            // Get rid of variable label incorporations
            string sequenceFixed = "";
            for (int i = 0; i < sequenceOriginal.Length; i++)
            {
                if (sequenceOriginal[i].Equals('l')) sequenceFixed += 'L';
                else if (sequenceOriginal[i].Equals('r')) sequenceFixed += 'R';
                else if (sequenceOriginal[i].Equals('k'))
                {
                    bool lysineSiteFound = false;
                    int count = 0;
                    while (!lysineSiteFound && count < lysineSites.Length)
                    {
                        for (int n = 0; n < lysineSites.Length; n++)
                        {
                            if (lysineSites[n] == i + 1 && lysinePTM[n]) lysineSiteFound = true;
                            count++;
                        }
                    }
                    if (!lysineSiteFound) sequenceFixed += 'K';
                    else sequenceFixed += 'k';
                }
                else sequenceFixed += sequenceOriginal[i];
            }
            
            // Initialize all peptide properties
            sequence = sequenceFixed;
            if (sequence.Contains('I'))
            {
                sequenceItoL = sequence.Replace('I', 'L');
            }
            else
            {
                sequenceItoL = sequence;
            }
            numChannels = Form1.NUMCHANNELS;
            numIsotopes = Form1.NUMISOTOPES;
            numIsotopologues = Form1.NUMISOTOPOLOGUES;
            numClusters = Form1.NUMCLUSTERS;

            // Deal with variable mods that do not affect quantification (e.g., oxidation, phosphorylation)
            sequenceNoMods = "";
            List<int> oxidationPositions = new List<int>();
            List<int> phosphorylationPositions = new List<int>();
            List<int> tyrosineNHSPositions = new List<int>();
            List<int> kAcPositions = new List<int>();
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
                    if (Form1.LABELSPERCHANNEL[0][0].LabelType.HasFlag(LabelType.NHS))
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
                    if (Form1.ACETYLATION)
                    {
                        kAcPositions.Add(i);
                    }
                    else if (Form1.UBIQUITIN)
                    {
                        kGGPositions.Add(i);
                    }
                }   
                else
                {
                    sequenceNoMods += sequence[i];
                }
            }

            peptide = new Peptide(sequenceNoMods);
            peptide.SetModification(NamedChemicalFormula.Carbamidomethyl, 'C');

            setTheoreticalIsotopicDistribution();

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
            //if (numIsotopologues > 1)
            //{
            //    firstSearchMassRange = new MassTolerance(MassToleranceType.PPM, 50.0);
            //    if (Form1.NHSISOTOPOLOGUE)
            //    {
            //        numIsotopologueLabels = countResidues('K', peptide.Sequence) + 1 + tyrosineNHSPositions.Count;
            //        firstSearchMassRange = new MassTolerance(MassToleranceType.PPM, 50.0);
            //    }
            //    else if (Form1.LEUISOTOPOLOGUE)
            //    {
            //        numIsotopologueLabels = countResidues('L', peptide.Sequence);
            //        firstSearchMassRange = new MassTolerance(MassToleranceType.PPM, 50.0);
            //    }
            //    else
            //    {
            //        numIsotopologueLabels = countResidues('K', peptide.Sequence);
            //        firstSearchMassRange = new MassTolerance(MassToleranceType.PPM, 50.0);

            //        if (numClusters > 1 || (numClusters == 1 && Form1.CLUSTERLABELS.Count > 0))
            //        {
            //            if (Form1.NHSCLUSTER)
            //            {
            //                numClusterLabels = countResidues('K', peptide.Sequence) + 1 + tyrosineNHSPositions.Count;
            //            }
            //        }
            //    }
            //}
            //else
            //{
            //    if (Form1.LEUCLUSTER)
            //    {
            //        numClusterLabels = countResidues('L', peptide.Sequence);
            //        firstSearchMassRange = new MassTolerance(MassToleranceType.PPM, 50.0);
            //    }
            //    else if (Form1.CYSCLUSTER)
            //    {
            //        numClusterLabels = countResidues('C', peptide.Sequence);
            //        firstSearchMassRange = new MassTolerance(MassToleranceType.PPM, 50.0);
            //    }
            //    else if (Form1.ARGCLUSTER)
            //    {
            //        numClusterLabels = countResidues('R', peptide.Sequence);
            //        firstSearchMassRange = new MassTolerance(MassToleranceType.PPM, 50.0);
            //    }
            //    else if (Form1.NHSCLUSTER)
            //    {
            //        numClusterLabels = countResidues('K', peptide.Sequence) + 1 + tyrosineNHSPositions.Count;
            //        firstSearchMassRange = new MassTolerance(MassToleranceType.PPM, 50.0);
            //    }
            //    else if (Form1.LYSCLUSTER)
            //    {
            //        numClusterLabels = countResidues('K', peptide.Sequence);
            //        firstSearchMassRange = new MassTolerance(MassToleranceType.PPM, 50.0);
            //    }
            //    else
            //    {
            //        numClusterLabels = 1;
            //        firstSearchMassRange = new MassTolerance(MassToleranceType.PPM, 50.0);
            //    }
            //}

            //theoMasses = new double[numChannels, 1];

            // Initialize all identification data members
            PSMs = new Dictionary<MSDataFile, List<PeptideSpectralMatch>>();
            PSMs.Add(rawFile, new List<PeptideSpectralMatch>());
            PSM = new PeptideSpectralMatch(scanNumber, rawFile, charge, eValue, fileNameID);
            PSMs.TryGetValue(rawFile, out psms);
            psms.Add(PSM);
            bestPSMs = new Dictionary<MSDataFile, PeptideSpectralMatch>();
            bestPSMs.Add(rawFile, PSM);
            firstSearchMassRange = new MassTolerance(MassToleranceType.PPM, 50.0);
            
            // Initialize all quantitation data members
            completeTotalIntensity = new double[numChannels, numIsotopes + 1];
            totalIntensity = new double[numChannels, numIsotopes + 1];
            if (numIsotopologues > 1)
            {
                heavyToLightRatioSum = new double[(numIsotopologues - 1) * numClusters];
                heavyToLightRatioMedian = new double[(numIsotopologues - 1) * numClusters, 1];
            }
            else
            {
                heavyToLightRatioSum = new double[(numChannels - 1)];
                heavyToLightRatioMedian = new double[numChannels - 1, 1];
            }
            allPairs = new Dictionary<MSDataFile, List<Pair>>();
            allPairs.Add(rawFile, new List<Pair>());
            completePairs = new Dictionary<MSDataFile, List<Pair>>();
            completePairs.Add(rawFile, new List<Pair>());
            missingChannelPattern = new int[numClusters, 2];
            numPeptideVersions = 1;

            peptideVersions = new Peptide[numChannels];
            theoMasses = new double[numChannels, 1];

            for (int i = 0; i < numChannels; i++)
            {
                peptideVersions[i] = new Peptide(peptide);

                if (phosphorylationPositions.Count > 0)
                {
                    foreach (int position in phosphorylationPositions)
                    {
                        peptideVersions[i].SetModification(NamedChemicalFormula.Phosphorylation, position + 1);
                    }
                }
                foreach (Label label in Form1.LABELSPERCHANNEL[i])
                {
                    peptideVersions[i].SetModification(label.LabelModification, label.LabelSites);

                    if (tyrosineNHSPositions.Count > 0 && label.LabelType.HasFlag(LabelType.NHS) && label.LabelSites.HasFlag(ModificationSites.NPep))
                    {
                        foreach (int position in tyrosineNHSPositions)
                        {
                            peptideVersions[i].SetModification(label.LabelModification, position + 1);
                        }
                    }

                    if (kAcPositions.Count > 0 && label.LabelType.HasFlag(LabelType.AminoAcid) && label.LabelSites.HasFlag(ModificationSites.K))
                    {
                        foreach (int position in kAcPositions)
                        {
                            peptideVersions[i].SetModification(new Mass(label.LabelModification.MonoisotopicMass + NamedChemicalFormula.Acetyl.MonoisotopicMass), position + 1);
                        }
                    }

                    if (kGGPositions.Count > 0 && label.LabelType.HasFlag(LabelType.AminoAcid) && label.LabelSites.HasFlag(ModificationSites.K))
                    {
                        foreach (int position in kGGPositions)
                        {
                            peptideVersions[i].SetModification(new Mass(label.LabelModification.MonoisotopicMass + 114.0429), position + 1);
                        }
                    }
                }
                theoMasses[i, 0] = peptideVersions[i].MonoisotopicMass;
                if (i > 0 && theoMasses[i,0] > theoMasses[i-1,0]) numPeptideVersions++;
            }

            if (numIsotopologues < 2)
            {
                isotopologueLabel = false;
                if (numPeptideVersions > 1)
                {
                    clusterLabel = true;
                    labeled = true;
                }
                else
                {
                    clusterLabel = false;
                    labeled = false;
                }
                
            }
            else
            {
                if (numPeptideVersions > 1)
                {
                    if (numClusters > 1)
                    {
                        if (numPeptideVersions == numChannels)
                        {
                            isotopologueLabel = true;
                            clusterLabel = true;
                            labeled = true;
                        }
                        else if (numPeptideVersions == numClusters)
                        {
                            isotopologueLabel = false;
                            clusterLabel = true;
                            labeled = true;
                        }
                        else if (numPeptideVersions == numIsotopologues)
                        {
                            isotopologueLabel = true;
                            clusterLabel = false;
                            labeled = true;
                        }
                    }
                    else
                    {
                        clusterLabel = false;
                        isotopologueLabel = true;
                        labeled = true;

                    }
                }
                else
                {
                    isotopologueLabel = false;
                    clusterLabel = false;
                    labeled = false;
                }
            }
            
            //// Set theoretical masses of each channel
            //if (numIsotopologueLabels == 0 && numIsotopologues > 1)
            //{
            //    numPeptideVersions = 1;
            //    peptideVersions = new Peptide[numPeptideVersions];
            //    peptideVersions[0] = new Peptide(peptide);
            //    if (Form1.CLUSTERLABELS.Count > 0)
            //    {
            //        peptideVersions[0].SetModification(Form1.CLUSTERLABELS[0], ModificationSites.NPep);   
            //    }
            //}
            //else if (numIsotopologues > 1 && numIsotopologueLabels > 0)
            //{
            //    if (numClusters > 1)
            //    {
            //        if (numClusterLabels > 0) numPeptideVersions = numIsotopologues * numClusters;
            //        else numPeptideVersions = numIsotopologues;
            //    }
            //    else
            //    {
            //        numPeptideVersions = numIsotopologues;
            //    }
            //    peptideVersions = new Peptide[numPeptideVersions];
                
            //    for (int i = 0; i < numPeptideVersions; i++)
            //    {
            //        peptideVersions[i] = new Peptide(peptide);
            //    }

            //    if ((numClusters == 1 && Form1.CLUSTERLABELS.Count == 0) || numClusterLabels == 0)
            //    {
            //        for (int i = 0; i < Form1.ISOTOPOLOGUELABELS.Count; i++)
            //        {
            //            if (Form1.LYSISOTOPOLOGUE)
            //            {
            //                peptideVersions[i].SetModification(Form1.ISOTOPOLOGUELABELS[i], ModificationSites.K);

            //                if (kGGPositions.Count > 0)
            //                {
            //                    foreach (int position in kGGPositions)
            //                    {
            //                        peptideVersions[i].SetModification(new Mass(Form1.ISOTOPOLOGUELABELS[i].MonoisotopicMass + 114.0429), position + 1);
            //                    }
            //                }
            //                else if (kAcPositions.Count > 0)
            //                {
            //                    foreach (int position in kAcPositions)
            //                    {
            //                        peptideVersions[i].SetModification(new Mass(Form1.ISOTOPOLOGUELABELS[i].MonoisotopicMass + 42.0106), position + 1);
            //                    }
            //                }
            //            }
            //            else if (Form1.LEUISOTOPOLOGUE)
            //            {
            //                peptideVersions[i].SetModification(Form1.ISOTOPOLOGUELABELS[i], ModificationSites.L);

            //                if (kGGPositions.Count > 0)
            //                {
            //                    foreach (int position in kGGPositions)
            //                    {
            //                        peptideVersions[i].SetModification(new Mass(114.0429), position + 1);
            //                    }
            //                }
            //                else if (kAcPositions.Count > 0)
            //                {
            //                    foreach (int position in kAcPositions)
            //                    {
            //                        peptideVersions[i].SetModification(new Mass(42.0106), position + 1);
            //                    }
            //                }
            //            }
            //            else if (Form1.NHSISOTOPOLOGUE)
            //            {
            //                peptideVersions[i].SetModification(Form1.ISOTOPOLOGUELABELS[i], ModificationSites.NPep | ModificationSites.K);

            //                if (tyrosineNHSPositions.Count > 0)
            //                {
            //                    foreach (int position in tyrosineNHSPositions)
            //                    {
            //                        peptideVersions[i].SetModification(Form1.ISOTOPOLOGUELABELS[i], position + 1);
            //                    }
            //                }
            //            }
            //        }
            //    }
            //    else if (numClusters == 1)
            //    {
            //        for (int i = 0; i < Form1.ISOTOPOLOGUELABELS.Count; i++)
            //        {
            //            if (Form1.LYSISOTOPOLOGUE)
            //            {
            //                peptideVersions[i].SetModification(Form1.ISOTOPOLOGUELABELS[i], ModificationSites.K);
            //                peptideVersions[i].SetModification(Form1.CLUSTERLABELS[0], ModificationSites.NPep);

            //                if (tyrosineNHSPositions.Count > 0)
            //                {
            //                    foreach (int position in tyrosineNHSPositions)
            //                    {
            //                        peptideVersions[i].SetModification(Form1.CLUSTERLABELS[0], position + 1);
            //                    }
            //                }
            //                else if (kGGPositions.Count > 0)
            //                {
            //                    foreach (int position in kGGPositions)
            //                    {
            //                        peptideVersions[i].SetModification(new Mass(Form1.ISOTOPOLOGUELABELS[i].MonoisotopicMass + 114.0429), position + 1);
            //                    }
            //                }
            //                else if (kAcPositions.Count > 0)
            //                {
            //                    foreach (int position in kAcPositions)
            //                    {
            //                        peptideVersions[i].SetModification(new Mass(Form1.ISOTOPOLOGUELABELS[i].MonoisotopicMass + 42.0106), position + 1);
            //                    }
            //                }
            //            }
            //        }
            //    }

            //    else if (numClusters > 2)
            //    {
            //        int channelIndex;

            //        for (int c = 0; c < numClusters; c++)
            //        {
            //            channelIndex = c * numIsotopologues;
            //            for (int i = channelIndex; i < channelIndex + numIsotopologues; i++)
            //            {
            //                if (Form1.ARGCLUSTER)
            //                {
            //                    if (c > 0)
            //                    {
            //                        peptideVersions[i].SetModification(Form1.CLUSTERLABELS[c - 1], ModificationSites.R);
            //                    }

            //                    if (Form1.LYSISOTOPOLOGUE)
            //                    {
            //                        peptideVersions[i].SetModification(Form1.ISOTOPOLOGUELABELS[i], ModificationSites.K);

            //                        if (kGGPositions.Count > 0)
            //                        {
            //                            foreach (int position in kGGPositions)
            //                            {
            //                                peptideVersions[i].SetModification(new Mass(Form1.ISOTOPOLOGUELABELS[i].MonoisotopicMass + 114.0429), position + 1);
            //                            }
            //                        }
            //                        else if (kAcPositions.Count > 0)
            //                        {
            //                            foreach (int position in kAcPositions)
            //                            {
            //                                peptideVersions[i].SetModification(new Mass(Form1.ISOTOPOLOGUELABELS[i].MonoisotopicMass + 42.0106), position + 1);
            //                            }
            //                        }
            //                    }
            //                    else if (Form1.LEUISOTOPOLOGUE)
            //                    {
            //                        peptideVersions[i].SetModification(Form1.ISOTOPOLOGUELABELS[i], ModificationSites.L);

            //                        if (kGGPositions.Count > 0)
            //                        {
            //                            foreach (int position in kGGPositions)
            //                            {
            //                                peptideVersions[i].SetModification(new Mass(114.0429), position + 1);
            //                            }
            //                        }
            //                        else if (kAcPositions.Count > 0)
            //                        {
            //                            foreach (int position in kAcPositions)
            //                            {
            //                                peptideVersions[i].SetModification(new Mass(42.0106), position + 1);
            //                            }
            //                        }
            //                    }
            //                }
            //                else if (Form1.LEUCLUSTER)
            //                {
            //                    if (c > 0)
            //                    {
            //                        peptideVersions[i].SetModification(Form1.CLUSTERLABELS[c - 1], ModificationSites.L);
            //                    }

            //                    if (Form1.LYSISOTOPOLOGUE)
            //                    {
            //                        peptideVersions[i].SetModification(Form1.ISOTOPOLOGUELABELS[i], ModificationSites.K);

            //                        if (kGGPositions.Count > 0)
            //                        {
            //                            foreach (int position in kGGPositions)
            //                            {
            //                                peptideVersions[i].SetModification(new Mass(Form1.ISOTOPOLOGUELABELS[i].MonoisotopicMass + 114.0429), position + 1);
            //                            }
            //                        }
            //                        else if (kAcPositions.Count > 0)
            //                        {
            //                            foreach (int position in kAcPositions)
            //                            {
            //                                peptideVersions[i].SetModification(new Mass(Form1.ISOTOPOLOGUELABELS[i].MonoisotopicMass + 42.0106), position + 1);
            //                            }
            //                        }
            //                    }
            //                }
            //                else if (Form1.NHSCLUSTER)
            //                {
            //                    peptideVersions[i].SetModification(Form1.CLUSTERLABELS[c], ModificationSites.NPep);
            //                    peptideVersions[i].SetModification(Form1.ISOTOPOLOGUELABELS[i], ModificationSites.K);

            //                    if (tyrosineNHSPositions.Count > 0)
            //                    {
            //                        foreach (int position in tyrosineNHSPositions)
            //                        {
            //                            peptideVersions[i].SetModification(Form1.CLUSTERLABELS[c], position);
            //                        }
            //                    }
            //                }
            //            }
            //        }
            //    }
            //}
            //else if (numIsotopologues < 2 && numClusterLabels > 0)
            //{
            //    numPeptideVersions = numClusters;
            //    peptideVersions = new Peptide[numPeptideVersions];

            //    for (int i = 0; i < numPeptideVersions; i++)
            //    {
            //        peptideVersions[i] = new Peptide(peptide);
            //    }

            //    for (int i = 0; i < Form1.CLUSTERLABELS.Count; i++)
            //    {
            //        if (Form1.LYSCLUSTER)
            //        {
            //            peptideVersions[i].SetModification(Form1.CLUSTERLABELS[i], ModificationSites.K);
            //        }
            //        else if (Form1.ARGCLUSTER)
            //        {
            //            peptideVersions[i + 1].SetModification(Form1.CLUSTERLABELS[i], ModificationSites.R);
            //        }
            //        else if (Form1.LEUCLUSTER)
            //        {
            //            peptideVersions[i + 1].SetModification(Form1.CLUSTERLABELS[i], ModificationSites.L);
            //        }
            //        else if (Form1.CYSCLUSTER)
            //        {
            //            peptideVersions[i + 1].SetModification(Form1.CLUSTERLABELS[i], ModificationSites.L);
            //        }
            //        else if (Form1.NHSCLUSTER)
            //        {
            //            peptideVersions[i].SetModification(Form1.CLUSTERLABELS[i], ModificationSites.NPep | ModificationSites.K);

            //            if (tyrosineNHSPositions.Count > 0)
            //            {
            //                foreach (int position in tyrosineNHSPositions)
            //                {
            //                    peptideVersions[i].SetModification(Form1.CLUSTERLABELS[i], position);
            //                }
            //            }

            //        }
            //        else
            //        {
            //            peptideVersions[i + 1].SetModification(Form1.CLUSTERLABELS[i], peptideVersions[i + 1].Length - 1);
            //        }
            //    }
            //}
            //else
            //{
            //    peptideVersions = new Peptide[numPeptideVersions];
            //    peptideVersions[0] = new Peptide(peptide);
            //}

            //theoMasses = new double[numPeptideVersions, 1];

            //if (peptideVersions != null)
            //{
            //    for (int i = 0; i < peptideVersions.Length; i++)
            //    {
            //        theoMasses[i, 0] = peptideVersions[i].MonoisotopicMass;
            //    }
            //}
                
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
            PeptideSpectralMatch best = null;
            foreach (KeyValuePair<MSDataFile, List<PeptideSpectralMatch>> psmSet in PSMs)
            {
                List<PeptideSpectralMatch> psms = psmSet.Value;
                List<PeptideSpectralMatch> eValueSortedPSMs = null;
                List<PeptideSpectralMatch> resolvabilitySortedPSMs = null;
                List<PeptideSpectralMatch> bestResolvableChargePSMs = null;
                List<PeptideSpectralMatch> eValueResolvabilitySortedPSMs = null;
                PeptideSpectralMatch currentPSM;
                double currentTheoResolvability;
                double currentExpResolvability;
                double currentResolvability;

                if (psms.Count == 1)
                {
                    bestPSMs[psmSet.Key] = psms[0];
                    if (best == null || bestPSMs[psmSet.Key].Resolvability > best.Resolvability || bestPSMs[psmSet.Key].EValue < best.EValue) best = bestPSMs[psmSet.Key];
                }
                else if ((psms.Count > 1 && numIsotopologues < 2) || (psms.Count > 1 && numIsotopologues > 1 && !isotopologueLabel))
                {
                    eValueSortedPSMs = psms.OrderBy(psm => psm.EValue).ToList();
                    bestPSMs[psmSet.Key] = eValueSortedPSMs[0];
                    if (best == null || bestPSMs[psmSet.Key].EValue < best.EValue) best = bestPSMs[psmSet.Key];
                }
                else if (psms.Count > 1 && isotopologueLabel)
                {
                    double coefficient = (Math.Sqrt(2 * Math.Log(100.0 / quantSeparation))) / (Math.Sqrt(2 * Math.Log(2)));
                    double spacing = spacingMassRange[1, 0].Mean;

                    for (int i = 0; i < psms.Count; i++)
                    {
                        currentPSM = psms[i];
                        currentTheoResolvability = coefficient * ((Mass.MzFromMass(theoMasses[0, 0], currentPSM.Charge) / (quantResolution * Math.Sqrt(400 / Mass.MzFromMass(theoMasses[0, 0], currentPSM.Charge)))));
                        currentExpResolvability = spacing / (double)currentPSM.Charge;
                        currentResolvability = currentExpResolvability / currentTheoResolvability;
                        currentPSM.Resolvability = currentResolvability;
                    }   
                 
                    resolvabilitySortedPSMs = psms.OrderByDescending(psm => psm.Resolvability).ToList();
                    int resolvableCharge = resolvabilitySortedPSMs[0].Charge;
                    bestResolvableChargePSMs = new List<PeptideSpectralMatch>();

                    for (int i = 0; i < resolvabilitySortedPSMs.Count; i++)
                    {
                        if (resolvabilitySortedPSMs[0].Charge == resolvableCharge) bestResolvableChargePSMs.Add(resolvabilitySortedPSMs[i]);
                    }

                    if (bestResolvableChargePSMs.Count > 1)
                    {
                        eValueResolvabilitySortedPSMs = bestResolvableChargePSMs.OrderBy(psm => psm.EValue).ToList();
                        bestPSMs[psmSet.Key] = eValueResolvabilitySortedPSMs[0];
                    }
                    else
                    {
                        bestPSMs[psmSet.Key] = bestResolvableChargePSMs[0];
                    }

                    if (best == null || bestPSMs[psmSet.Key].Resolvability > best.Resolvability || bestPSMs[psmSet.Key].EValue < best.EValue) best = bestPSMs[psmSet.Key];
                }
            }
            bestPSM = best;
        }
        
        /* Calculates the MS1 scan range to use for quantification
         * Includes all PSMs, extending above and below according to the entered retention time window
         */
        public void calculateScanRange(MSDataFile rawFile, double rTWindowMin, double rTWindowMax)
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
                rTRange = new Range<double>(firstScanTime - rTWindowMin, lastScanTime + rTWindowMax);
            }
            else
            {
                rTRange = new Range<double>(bestScanTime - rTWindowMin, bestScanTime + rTWindowMax);
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
                    if (newPeak.GetSignalToNoise() < minSignalToNoise || newPeak.Y > maxRawIntensity)
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

                    if (mZSortedPeaks.Count == 2)
                    {
                        Spacings.Add(new Spacing(sequence, CurrentScan.RetentionTime, CurrentScan.SpectrumNumber, Charge, countResidues('K', sequence), Isotope, this, mZSortedPeaks[0], mZSortedPeaks[1]));
                    }                    
                }
                convertedPeaks = cleanPeaks(convertedPeaks, Charge);
            }
            return convertedPeaks;
        }

        /* Checks to see if search window should be shifted based on the detected peaks and the errors
         * 
         */
        public MassRange checkMappedPeaks(ILabeledPeak[] MappedPeaks, double ToleranceWidth, MassTolerance Tolerance, int Cluster, int Isotope, int Charge)
        {
            MassRange newSearchWindow = null;
            List<int> peakPositions = new List<int>();
            List<double> peakPpmErrors = new List<double>();

            for (int i = 0; i < MappedPeaks.Length; i++)
            {
                int channelIndex = (Cluster * numIsotopologues) + i;
                if (MappedPeaks[i] != null && MappedPeaks[i].GetSignalToNoise() >= Form1.MINIMUMSN)
                {
                    peakPositions.Add(i);
                    peakPpmErrors.Add(Math.Abs(MassTolerance.GetTolerance(Mass.MassFromMz(MappedPeaks[i].X, Charge), adjustedTheoMasses[channelIndex, Isotope], MassToleranceType.PPM)));
                }
            }

            int peakCount = peakPositions.Count;

            if (peakCount >= peaksNeeded)
            {
                double newMinMZ = 0;
                double newMaxMZ = 0;
                if (numIsotopologues == 2)
                {
                    // Lightest peak detected
                    if (peakPositions[0] == 0)
                    {
                        // Ppm error is outside range
                        if (peakPpmErrors[0] >= ppmTolerance.Value)
                        {
                            newMinMZ = new MassRange(MappedPeaks[0].X, Tolerance).Minimum;
                            newMaxMZ = newMinMZ + ToleranceWidth;
                        }
                        // Shift search window to lower m/z
                        else
                        {
                            newMaxMZ = new MassRange(MappedPeaks[0].X, Tolerance).Maximum;
                            newMinMZ = newMaxMZ - ToleranceWidth;
                        }
                    }
                    // Heaviest peak detected
                    else if (peakPositions[0] == 1)
                    {
                        // Ppm error is outside range
                        if (peakPpmErrors[0] >= ppmTolerance.Value)
                        {
                            newMaxMZ = new MassRange(MappedPeaks[1].X, Tolerance).Maximum;
                            newMinMZ = newMaxMZ - ToleranceWidth;
                        }
                        // Shift search window to higher m/z
                        else
                        {
                            newMinMZ = new MassRange(MappedPeaks[1].X, Tolerance).Minimum;
                            newMaxMZ = newMinMZ + ToleranceWidth;
                        }
                    }
                }
                else if (numIsotopologues == 3)
                {
                    if (peakCount == 1)
                    {
                        // Lightest peak detected
                        if (peakPositions[0] == 0)
                        {
                            if (peakPpmErrors[0] >= ppmTolerance.Value)
                            {
                                newMinMZ = new MassRange(MappedPeaks[0].X, Tolerance).Minimum;
                                newMaxMZ = newMinMZ + ToleranceWidth;
                            }
                            else
                            {
                                newMaxMZ = new MassRange(MappedPeaks[0].X, Tolerance).Maximum;
                                newMinMZ = newMaxMZ - ToleranceWidth;
                            }
                        }
                        // Heaviest peak detected
                        else if (peakPositions[0] == 2)
                        {
                            if (peakPpmErrors[0] >= ppmTolerance.Value)
                            {
                                newMaxMZ = new MassRange(MappedPeaks[2].X, Tolerance).Maximum;
                                newMinMZ = newMaxMZ - ToleranceWidth;
                            }
                            else
                            {
                                newMinMZ = new MassRange(MappedPeaks[2].X, Tolerance).Minimum;
                                newMaxMZ = newMinMZ + ToleranceWidth;
                            }
                        }
                    }
                    else if (peakCount == 2)
                    {
                        // Lightest two peaks detected
                        if (peakPositions[0] == 0 && peakPositions[1] == 1)
                        {
                            if (peakPpmErrors[0] >= ppmTolerance.Value || peakPpmErrors[1] >= ppmTolerance.Value)
                            {
                                newMinMZ = new MassRange(MappedPeaks[0].X, Tolerance).Minimum;
                                newMaxMZ = newMinMZ + ToleranceWidth;
                            }
                            else
                            {
                                newMaxMZ = new MassRange(MappedPeaks[1].X, Tolerance).Maximum;
                                newMinMZ = newMaxMZ - ToleranceWidth;
                            }
                        }
                        // Heaviest two peaks detected
                        else if (peakPositions[0] == 1 && peakPositions[1] == 2)
                        {
                            if (peakPpmErrors[0] >= ppmTolerance.Value || peakPpmErrors[1] >= ppmTolerance.Value)
                            {
                                newMaxMZ = new MassRange(MappedPeaks[2].X, Tolerance).Maximum;
                                newMinMZ = newMaxMZ - ToleranceWidth;
                            }
                            else
                            {
                                newMinMZ = new MassRange(MappedPeaks[1].X, Tolerance).Minimum;
                                newMaxMZ = newMinMZ + ToleranceWidth;
                            }
                        }
                    }
                }
                else if (numIsotopologues == 4)
                {
                    if (peakCount == 2)
                    {
                        // Lightest peaks detected
                        if (peakPositions[0] == 0 && peakPositions[1] == 1)
                        {
                            if (peakPpmErrors[0] >= ppmTolerance.Value || peakPpmErrors[1] >= ppmTolerance.Value)
                            {
                                newMinMZ = new MassRange(MappedPeaks[0].X, Tolerance).Minimum;
                                newMaxMZ = newMinMZ + ToleranceWidth;
                            }
                            else
                            {
                                newMaxMZ = new MassRange(MappedPeaks[1].X, Tolerance).Maximum;
                                newMinMZ = newMaxMZ - ToleranceWidth;
                            }
                        }
                        else if(peakPositions[0] == 0 && peakPositions[1] == 2)
                        {
                            if (peakPpmErrors[0] >= ppmTolerance.Value || peakPpmErrors[1] >= ppmTolerance.Value)
                            {
                                newMinMZ = new MassRange(MappedPeaks[0].X, Tolerance).Minimum;
                                newMaxMZ = newMinMZ + ToleranceWidth;
                            }
                            else
                            {
                                newMaxMZ = new MassRange(MappedPeaks[2].X, Tolerance).Maximum;
                                newMinMZ = newMaxMZ - ToleranceWidth;
                            }
                        }
                        // Heaviest peaks detected
                        else if (peakPositions[0] == 2 && peakPositions[1] == 3)
                        {
                            if (peakPpmErrors[0] >= ppmTolerance.Value || peakPpmErrors[1] >= ppmTolerance.Value)
                            {
                                newMaxMZ = new MassRange(MappedPeaks[3].X, Tolerance).Maximum;
                                newMinMZ = newMaxMZ - ToleranceWidth;
                            }
                            else
                            {
                                newMinMZ = new MassRange(MappedPeaks[2].X, Tolerance).Minimum;
                                newMaxMZ = newMinMZ + ToleranceWidth;
                            }
                        }
                        else if (peakPositions[0] == 1 && peakPositions[1] == 3)
                        {
                            if (peakPpmErrors[0] >= ppmTolerance.Value || peakPpmErrors[1] >= ppmTolerance.Value)
                            {
                                newMaxMZ = new MassRange(MappedPeaks[3].X, Tolerance).Maximum;
                                newMinMZ = newMaxMZ - ToleranceWidth;
                            }
                            else
                            {
                                newMinMZ = new MassRange(MappedPeaks[1].X, Tolerance).Minimum;
                                newMaxMZ = newMinMZ + ToleranceWidth;
                            }
                        }
                    }
                    else if (peakCount == 3)
                    {
                        if (peakPositions[0] == 0 && peakPositions[1] == 1 && peakPositions[2] == 2)
                        {
                            if (peakPpmErrors[0] >= ppmTolerance.Value || peakPpmErrors[1] >= ppmTolerance.Value || peakPpmErrors[2] >= ppmTolerance.Value)
                            {
                                newMinMZ = new MassRange(MappedPeaks[0].X, Tolerance).Minimum;
                                newMaxMZ = newMinMZ + ToleranceWidth;
                            }
                            else
                            {
                                newMaxMZ = new MassRange(MappedPeaks[2].X, Tolerance).Maximum;
                                newMinMZ = newMaxMZ - ToleranceWidth;
                            }
                        }
                        else if (peakPositions[0] == 1 && peakPositions[1] == 2 && peakPositions[2] == 3)
                        {
                            if (peakPpmErrors[0] >= ppmTolerance.Value || peakPpmErrors[1] >= ppmTolerance.Value || peakPpmErrors[2] >= ppmTolerance.Value)
                            {
                                newMaxMZ = new MassRange(MappedPeaks[2].X, Tolerance).Maximum;
                                newMinMZ = newMaxMZ - ToleranceWidth;
                            }
                            else
                            {
                                newMinMZ = new MassRange(MappedPeaks[1].X, Tolerance).Minimum;
                                newMaxMZ = newMinMZ + ToleranceWidth;
                            }
                        }
                    }
                }
                else if (numIsotopologues == 6)
                {

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
            ILabeledPeak[] peaksFound = new ILabeledPeak[numChannels];
            double injectionTime = rawFile.GetInjectionTime(current.SpectrumNumber);
            Pair pair = new Pair(this, rawFile, current.SpectrumNumber, charge, injectionTime, current.RetentionTime);
            int peaksCount = 0;

            if (numIsotopologues > 1 && !Form1.CROSSCLUSTERQUANT)
            {
                for (int c = 0; c < numClusters; c++) // Start cluster loop
                {
                    if (GetTheoreticalResolvability(c))
                    {
                        for (int j = 0; j < numIsotopes; j++) // Start isotope loop
                        {
                            peaksFound = new ILabeledPeak[numChannels];
                            peaksCount = 0;
                            int channelIndex = c * numIsotopologues;

                            // Start pattern search
                            MassTolerance tolerance = new MassTolerance(MassToleranceType.PPM, ppmTolerance.Value);

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
                                        peaksFound[channelIndex + i] = orderedByMzPeaks[i];
                                    }
                                }
                                else
                                {
                                    //if (!Form1.NOISEBANDCAP)
                                    //{
                                    //    if (topSignalToNoisePeaks.Count >= numIsotopologues / 2)
                                    //    {
                                    //        ILabeledPeak[] mappedPeaks = mapPeaks(topSignalToNoisePeaks, channelIndex, channelIndex + (numIsotopologues - 1), j, charge);
                                    //        MassRange newSearchWindow = checkMappedPeaks(mappedPeaks, totalRange.Width, tolerance, c, j, charge);

                                    //        if (newSearchWindow != null)
                                    //        {
                                    //            List<ILabeledPeak> secondTryPeaks = findPeakPatterns(newSearchWindow, charge, current);
                                    //            List<ILabeledPeak> secondTryTopSignalToNoisePeaks = null;
                                    //            List<ILabeledPeak> secondTryOrderedByMzPeaks = null;
                                    //            if (secondTryPeaks != null && secondTryPeaks.Count == numIsotopologues)
                                    //            {
                                    //                secondTryTopSignalToNoisePeaks = secondTryPeaks.OrderByDescending(peakSort => peakSort.GetSignalToNoise()).Take(numIsotopologues).ToList();
                                    //                secondTryOrderedByMzPeaks = secondTryTopSignalToNoisePeaks.OrderBy(peakSort => peakSort.X).ToList();
                                    //                for (int i = 0; i < numIsotopologues; i++)
                                    //                {
                                    //                    peaksFound[channelIndex + i, j] = secondTryOrderedByMzPeaks[i];
                                    //                }
                                    //            }
                                    //            else
                                    //            {
                                    //                for (int m = channelIndex; m < channelIndex + numIsotopologues; m++)
                                    //                {
                                    //                    peaksFound[m, j] = mappedPeaks[m - channelIndex];
                                    //                }
                                    //            }
                                    //        }
                                    //        else
                                    //        {
                                    //            for (int i = 0; i < numIsotopologues; i++)
                                    //            {
                                    //                peaksFound[channelIndex + i, j] = null;
                                    //            }
                                    //        }
                                    //    }
                                    //    else
                                    //    {
                                    //        for (int i = 0; i < numIsotopologues; i++)
                                    //        {
                                    //            peaksFound[channelIndex + i, j] = null;
                                    //        }
                                    //    }
                                    //}
                                    if (Form1.NOISEBANDCAP)
                                    {
                                        // Peaks found for an incomplete set of isotopologues by pattern detection
                                        if (topSignalToNoisePeaks.Count >= peaksNeeded)
                                        {
                                            ILabeledPeak[] mappedPeaks = mapPeaks(topSignalToNoisePeaks, channelIndex, channelIndex + (numIsotopologues - 1), j, charge);
                                            for (int m = channelIndex; m < channelIndex + numIsotopologues; m++)
                                            {
                                                peaksFound[m] = mappedPeaks[m - channelIndex];
                                            }
                                            //MassRange newSearchWindow = checkMappedPeaks(mappedPeaks, totalRange.Width, tolerance, c, j, charge);

                                            //if (newSearchWindow == null)
                                            //{
                                            //    for (int m = channelIndex; m < channelIndex + numIsotopologues; m++)
                                            //    {
                                            //        peaksFound[m, j] = mappedPeaks[m - channelIndex];
                                            //    }
                                            //}
                                            //else
                                            //{
                                            //    List<ILabeledPeak> secondTryPeaks = findPeakPatterns(newSearchWindow, charge, current);
                                            //    List<ILabeledPeak> secondTryTopSignalToNoisePeaks = null;
                                            //    List<ILabeledPeak> secondTryOrderedByMzPeaks = null;
                                            //    if (secondTryPeaks != null && secondTryPeaks.Count > topSignalToNoisePeaks.Count)
                                            //    {
                                            //        secondTryTopSignalToNoisePeaks = secondTryPeaks.OrderByDescending(peakSort => peakSort.GetSignalToNoise()).Take(numIsotopologues).ToList();
                                            //        if (secondTryTopSignalToNoisePeaks.Count == numIsotopologues)
                                            //        {
                                            //            secondTryOrderedByMzPeaks = secondTryTopSignalToNoisePeaks.OrderBy(peakSort => peakSort.X).ToList();
                                            //            for (int i = 0; i < numIsotopologues; i++)
                                            //            {
                                            //                peaksFound[channelIndex + i, j] = secondTryOrderedByMzPeaks[i];
                                            //            }
                                            //        }
                                            //        else if (secondTryTopSignalToNoisePeaks.Count >= peaksNeeded)
                                            //        {
                                            //            double additionalPPMError = MassTolerance.GetTolerance(newSearchWindow.Minimum, totalRange.Minimum, MassToleranceType.PPM);
                                            //            ILabeledPeak[] secondTryMappedPeaks = mapPeaks(secondTryTopSignalToNoisePeaks, channelIndex, channelIndex + (numIsotopologues - 1), j, charge, additionalPPMError);
                                            //            for (int m = channelIndex; m < channelIndex + numIsotopologues; m++)
                                            //            {
                                            //                peaksFound[m, j] = secondTryMappedPeaks[m - channelIndex];
                                            //            }
                                            //        }
                                            //    }
                                            //    else
                                            //    {
                                            //        for (int m = channelIndex; m < channelIndex + numIsotopologues; m++)
                                            //        {
                                            //            peaksFound[m, j] = mappedPeaks[m - channelIndex];
                                            //        }
                                            //    }
                                            //}
                                        }
                                    }
                                }
                            }

                            // Count detected peaks
                            for (int i = channelIndex; i < channelIndex + numIsotopologues; i++)
                            {
                                if (peaksFound[i] != null)
                                {
                                    peaksCount++;
                                }
                            }

                            // If enough peaks are detected, increase the cluster count
                            if (peaksCount >= peaksNeeded)
                            {
                                IsotopePair isotopePair = new IsotopePair(pair, j);    
                                for (int i = channelIndex; i < channelIndex + numIsotopologues; i++)
                                {
                                    isotopePair.ChannelPeaks[i] = peaksFound[i];
                                }
                                pair.IsotopePairs.Add(isotopePair);
                            }
                            // If not, set all peaks to null
                            else
                            {
                                for (int i = channelIndex; i < channelIndex + numIsotopologues; i++)
                                {
                                    peaksFound[i] = null;
                                }
                            }
                        } // End isotope loop
                    }
                } // End cluster loop

                //if (Form1.SEGMENTEDINJECTIONTIMES)
                //{
                //    bool peakFound = false;
                //    int channelCount = 0;
                //    double mZ = 0;
                //    double timeSum = 0;

                //    while (channelCount < numChannels && !peakFound)
                //    {
                //        for (int j = 0; j < numIsotopes; j++)
                //        {
                //            if (peaksFound[channelCount, j] != null)
                //            {
                //                mZ = peaksFound[channelCount, j].X;
                //                peakFound = true;
                //            }
                //        }                    
                //        channelCount++;
                //    }
                        
                //    Dictionary<int, Dictionary<Range<double>, double>> scanNumbers = null;
                //    if (Form1.INJECTIONTIMES.TryGetValue(rawFile.Name, out scanNumbers))
                //    {
                //        Dictionary<Range<double>, double> injectionTimes = null;

                //        if (scanNumbers.TryGetValue(current.SpectrumNumber, out injectionTimes))
                //        {
                //            List<double> times = new List<double>();
                //            foreach (KeyValuePair<Range<double>, double> segment in injectionTimes)
                //            {
                //                if (segment.Key.Contains(mZ))
                //                {
                //                    times.Add(segment.Value);
                //                    timeSum += segment.Value;
                //                }
                //            }

                //            pair.InjectionTime = (timeSum / ((double)times.Count)) / 1000.0;
                //        }
                //    }
                //}
                if (pair.IsotopePairs.Count > 0)
                {
                    allPairs[rawFile].Add(pair);
                }   
            }
            //else if (numIsotopologues > 1 && Form1.CROSSCLUSTERQUANT)
            //{
            //    if (allClustersResolvable)
            //    {
            //        int channelIndex;
            //        isotopeCount = 0;
            //        for (int j = 0; j < numIsotopes; j++) // Start cluster loop
            //        {
            //            peaksCount = 0;
            //            for (int c = 0; c < numClusters; c++) // Start isotope loop
            //            {
            //                channelIndex = c * numIsotopologues;

            //                // Start pattern search
            //                MassTolerance tolerance = new MassTolerance(MassToleranceType.PPM, ppmTolerance.Value);

            //                MassRange minRange = new MassRange(adjustedTheoMasses[channelIndex, j], tolerance);
            //                MassRange maxRange = new MassRange(adjustedTheoMasses[channelIndex + (numIsotopologues - 1), j], tolerance);

            //                double minMZ = Mass.MzFromMass(minRange.Minimum, charge);
            //                double maxMZ = Mass.MzFromMass(maxRange.Maximum, charge);

            //                MassRange totalRange = new MassRange(minMZ, maxMZ);

            //                List<ILabeledPeak> cleanedPeaks = findPeakPatterns(totalRange, charge, current, spacings, j);
            //                List<ILabeledPeak> topSignalToNoisePeaks = null;
            //                List<ILabeledPeak> orderedByMzPeaks = null;
            //                if (cleanedPeaks != null && cleanedPeaks.Count > 0)
            //                {
            //                    topSignalToNoisePeaks = cleanedPeaks.OrderByDescending(peakSort => peakSort.GetSignalToNoise()).Take(numIsotopologues).ToList();
            //                    // Peaks found for a complete set of isotopologues                             
            //                    if (topSignalToNoisePeaks.Count == numIsotopologues)
            //                    {
            //                        orderedByMzPeaks = topSignalToNoisePeaks.OrderBy(peakSort => peakSort.X).ToList();
            //                        for (int i = 0; i < numIsotopologues; i++)
            //                        {
            //                            peaksFound[channelIndex + i, j] = orderedByMzPeaks[i];
            //                        }
            //                    }
            //                    else
            //                    {
            //                        //if (!Form1.NOISEBANDCAP)
            //                        //{
            //                        //    if (topSignalToNoisePeaks.Count >= numIsotopologues / 2)
            //                        //    {
            //                        //        ILabeledPeak[] mappedPeaks = mapPeaks(topSignalToNoisePeaks, channelIndex, channelIndex + (numIsotopologues - 1), j, charge);
            //                        //        MassRange newSearchWindow = checkMappedPeaks(mappedPeaks, totalRange.Width, tolerance, c, j, charge);

            //                        //        if (newSearchWindow != null)
            //                        //        {
            //                        //            List<ILabeledPeak> secondTryPeaks = findPeakPatterns(newSearchWindow, charge, current);
            //                        //            List<ILabeledPeak> secondTryTopSignalToNoisePeaks = null;
            //                        //            List<ILabeledPeak> secondTryOrderedByMzPeaks = null;
            //                        //            if (secondTryPeaks != null && secondTryPeaks.Count == numIsotopologues)
            //                        //            {
            //                        //                secondTryTopSignalToNoisePeaks = secondTryPeaks.OrderByDescending(peakSort => peakSort.GetSignalToNoise()).Take(numIsotopologues).ToList();
            //                        //                secondTryOrderedByMzPeaks = secondTryTopSignalToNoisePeaks.OrderBy(peakSort => peakSort.X).ToList();
            //                        //                for (int i = 0; i < numIsotopologues; i++)
            //                        //                {
            //                        //                    peaksFound[channelIndex + i, j] = secondTryOrderedByMzPeaks[i];
            //                        //                }
            //                        //            }
            //                        //            else
            //                        //            {
            //                        //                for (int m = channelIndex; m < channelIndex + numIsotopologues; m++)
            //                        //                {
            //                        //                    peaksFound[m, j] = mappedPeaks[m - channelIndex];
            //                        //                }
            //                        //            }
            //                        //        }
            //                        //        else
            //                        //        {
            //                        //            for (int i = 0; i < numIsotopologues; i++)
            //                        //            {
            //                        //                peaksFound[channelIndex + i, j] = null;
            //                        //            }
            //                        //        }
            //                        //    }
            //                        //    else
            //                        //    {
            //                        //        for (int i = 0; i < numIsotopologues; i++)
            //                        //        {
            //                        //            peaksFound[channelIndex + i, j] = null;
            //                        //        }
            //                        //    }
            //                        //}
            //                        if (Form1.NOISEBANDCAP)
            //                        {
            //                            // Peaks found for an incomplete set of isotopologues by pattern detection
            //                            if (topSignalToNoisePeaks.Count >= 0)
            //                            {
            //                                ILabeledPeak[] mappedPeaks = mapPeaks(topSignalToNoisePeaks, channelIndex, channelIndex + (numIsotopologues - 1), j, charge);
            //                                MassRange newSearchWindow = checkMappedPeaks(mappedPeaks, totalRange.Width, tolerance, c, j, charge);

            //                                if (newSearchWindow == null)
            //                                {
            //                                    for (int m = channelIndex; m < channelIndex + numIsotopologues; m++)
            //                                    {
            //                                        peaksFound[m, j] = mappedPeaks[m - channelIndex];
            //                                    }
            //                                }
            //                                else
            //                                {
            //                                    List<ILabeledPeak> secondTryPeaks = findPeakPatterns(newSearchWindow, charge, current);
            //                                    List<ILabeledPeak> secondTryTopSignalToNoisePeaks = null;
            //                                    List<ILabeledPeak> secondTryOrderedByMzPeaks = null;
            //                                    if (secondTryPeaks != null && secondTryPeaks.Count > topSignalToNoisePeaks.Count)
            //                                    {
            //                                        secondTryTopSignalToNoisePeaks = secondTryPeaks.OrderByDescending(peakSort => peakSort.GetSignalToNoise()).Take(numIsotopologues).ToList();
            //                                        if (secondTryTopSignalToNoisePeaks.Count == numIsotopologues)
            //                                        {
            //                                            secondTryOrderedByMzPeaks = secondTryTopSignalToNoisePeaks.OrderBy(peakSort => peakSort.X).ToList();
            //                                            for (int i = 0; i < numIsotopologues; i++)
            //                                            {
            //                                                peaksFound[channelIndex + i, j] = secondTryOrderedByMzPeaks[i];
            //                                            }
            //                                        }
            //                                        else if (secondTryTopSignalToNoisePeaks.Count >= 0)
            //                                        {
            //                                            double additionalPPMError = MassTolerance.GetTolerance(newSearchWindow.Minimum, totalRange.Minimum, MassToleranceType.PPM);
            //                                            ILabeledPeak[] secondTryMappedPeaks = mapPeaks(secondTryTopSignalToNoisePeaks, channelIndex, channelIndex + (numIsotopologues - 1), j, charge, additionalPPMError);
            //                                            for (int m = channelIndex; m < channelIndex + numIsotopologues; m++)
            //                                            {
            //                                                peaksFound[m, j] = secondTryMappedPeaks[m - channelIndex];
            //                                            }
            //                                        }
            //                                    }
            //                                    else
            //                                    {
            //                                        for (int m = channelIndex; m < channelIndex + numIsotopologues; m++)
            //                                        {
            //                                            peaksFound[m, j] = mappedPeaks[m - channelIndex];
            //                                        }
            //                                    }
            //                                }
            //                            }
            //                        }
            //                    }

            //                    for (int m = channelIndex; m < channelIndex + numIsotopologues; m++)
            //                    {
            //                        if (peaksFound[m, j] != null) peaksCount++;
            //                    }
            //                }
            //            } // End cluster loop
            //            if (peaksCount >= peaksNeeded) isotopeCount++;
            //            else
            //            {
            //                for (int i = 0; i < numChannels; i++)
            //                {
            //                    peaksFound[i, j] = null;
            //                }
            //            }
            //        } // End isotope loop

            //        if (isotopeCount > 0)
            //        {
            //            pair = new Pair(this, rawFile, current.SpectrumNumber, injectionTime, current.RetentionTime);
            //            pair.peaks = peaksFound;
            //            if (pair.totalPeakCount > 0)
            //            {
            //                allPairs[rawFile].Add(pair);
            //            }
            //        }
            //    } 
            //}
            // Traditional SILAC (i.e., one channel per cluster)
            else
            {
                for (int j = 0; j < numIsotopes; j++)
                {
                    peaksCount = 0;
                    peaksFound = new ILabeledPeak[numChannels];
                    for (int i = 0; i < numChannels; i++)
                    {
                        peaksFound[i] = largestPeak(adjustedTheoMasses[i, j], current, ppmTolerance, rawFile);
                        if (peaksFound[i] != null)
                        {
                            peaksCount++;
                        }
                    }
                    if (peaksCount >= peaksNeeded)
                    {
                        IsotopePair isotopePair = new IsotopePair(pair, j);
                        for (int i = 0; i < numChannels; i++)
                        {
                            isotopePair.ChannelPeaks[i] = peaksFound[i];
                        }
                        pair.IsotopePairs.Add(isotopePair);
                        
                    }
                    else
                    {
                        for (int k = 0; k < numChannels; k++)
                        {
                            peaksFound[k] = null;
                        }
                    }
                }

                if (pair.IsotopePairs.Count > 0)
                {
                    allPairs[rawFile].Add(pair);
                }
            }
        }

        public List<ILabeledPeak> cleanPeaks(List<ILabeledPeak> initialList, int charge)
        {
            if (initialList != null && initialList.Count > 1)
            {
                double theoTh = (spacingMassRange[1, 0].Maximum * 0.9 * 1000.0) /  (double)charge;
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
                    error1 = Math.Abs(MassTolerance.GetTolerance(neutralMass, additionallyAdjustedTheoMasses[channelStart, isotope], MassToleranceType.PPM));
                    error2 = Math.Abs(MassTolerance.GetTolerance(neutralMass, additionallyAdjustedTheoMasses[channelEnd, isotope], MassToleranceType.PPM));
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
                    error1 = Math.Abs(MassTolerance.GetTolerance(neutralMass, additionallyAdjustedTheoMasses[channelStart, isotope], MassToleranceType.PPM));
                    error2 = Math.Abs(MassTolerance.GetTolerance(neutralMass, additionallyAdjustedTheoMasses[channelStart + 1, isotope], MassToleranceType.PPM));
                    error3 = Math.Abs(MassTolerance.GetTolerance(neutralMass, additionallyAdjustedTheoMasses[channelStart + 2, isotope], MassToleranceType.PPM));

                    // 3 possibilities
                    if (orderedByMzPeaks.Count == 1)
                    {
                        if (error1 < error2 && error1 < error3) mappedPeaks[0] = currentPeak;
                        else if (error2 < error1 && error2 < error3) mappedPeaks[1] = currentPeak;
                        else if (error3 < error1 && error3 < error2) mappedPeaks[2] = currentPeak;

                        return mappedPeaks;
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
                    error1 = Math.Abs(MassTolerance.GetTolerance(neutralMass, additionallyAdjustedTheoMasses[channelStart, isotope], MassToleranceType.PPM));
                    error2 = Math.Abs(MassTolerance.GetTolerance(neutralMass, additionallyAdjustedTheoMasses[channelStart + 1, isotope], MassToleranceType.PPM));
                    error3 = Math.Abs(MassTolerance.GetTolerance(neutralMass, additionallyAdjustedTheoMasses[channelStart + 2, isotope], MassToleranceType.PPM));
                    error4 = Math.Abs(MassTolerance.GetTolerance(neutralMass, additionallyAdjustedTheoMasses[channelStart + 3, isotope], MassToleranceType.PPM));

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
                    error1 = Math.Abs(MassTolerance.GetTolerance(neutralMass, additionallyAdjustedTheoMasses[channelStart, isotope], MassToleranceType.PPM));
                    error2 = Math.Abs(MassTolerance.GetTolerance(neutralMass, additionallyAdjustedTheoMasses[channelStart + 1, isotope], MassToleranceType.PPM));
                    error3 = Math.Abs(MassTolerance.GetTolerance(neutralMass, additionallyAdjustedTheoMasses[channelStart + 2, isotope], MassToleranceType.PPM));
                    error4 = Math.Abs(MassTolerance.GetTolerance(neutralMass, additionallyAdjustedTheoMasses[channelStart + 3, isotope], MassToleranceType.PPM));
                    error5 = Math.Abs(MassTolerance.GetTolerance(neutralMass, additionallyAdjustedTheoMasses[channelStart + 4, isotope], MassToleranceType.PPM));
                    error6 = Math.Abs(MassTolerance.GetTolerance(neutralMass, additionallyAdjustedTheoMasses[channelStart + 5, isotope], MassToleranceType.PPM));

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
        //public void checkPairCoalescence(MSDataFile rawFile)
        //{
        //    PeptideSpectralMatch best = PSMs[rawFile][0];
        //    int scanNumber = best.ScanNumber;
        //    double injectionTime;
        //    List<Pair> pairs;
        //    Pair pair;
        //    ILabeledPeak peak1;
        //    ILabeledPeak peak2;
        //    ILabeledPeak peak3;
        //    ILabeledPeak peak4;
        //    IsotopePair isotopePair;
        //    List<IsotopePair> lowerIntensity;
        //    List<IsotopePair> higherIntensity;
        //    double lowerFrequency;
        //    double higherFrequency;

        //    allPairs.TryGetValue(rawFile, out pairs);
        //    for (int p = 0; p < pairs.Count(); p++)
        //    {
        //        pair = pairs[p];
        //        injectionTime = rawFile.GetInjectionTime(pair.ScanNumber);
        //        for (int c = 0; c < numChannels; c += numIsotopologues)
        //        {
        //            for (int j = 0; j < numIsotopes; j++)
        //            {
        //                if (numIsotopologues == 2)
        //                {
        //                    peak1 = pair.peaks[c, j];
        //                    peak2 = pair.peaks[c + 1, j];

        //                    // Do not alter complete or null pairs or pairs with a maximum intensity below the set threshold
        //                    if (peak1 != null && peak2 != null)
        //                    {

        //                    }
        //                    else if (peak1 == null && peak2 == null)
        //                    {
                                
        //                    }
        //                    //else if (log10MaxIntensity <= System.Math.Log10(Form1.MAXIMUMDNL))
        //                    //{
                                
        //                    //}
        //                    // Proceed with the analysis for all other pairs
        //                    else
        //                    {
        //                        ILabeledPeak singlePeak;
        //                        ILabeledPeak otherPeak1;
        //                        ILabeledPeak otherPeak2;
        //                        double singlePeakIntensity;
        //                        Pair otherPair;

        //                        // Find the single peak
        //                        if (peak1 == null)
        //                        {
        //                            singlePeak = peak2;
        //                        }
        //                        else
        //                        {
        //                            singlePeak = peak1;
        //                        }
        //                        singlePeakIntensity = singlePeak.GetDenormalizedIntensity(injectionTime);

        //                        // Sort the peptide's other pairs (complete and incomplete) based on their intensities relative to the single peak
        //                        lowerIntensity = new List<IsotopePair>();
        //                        higherIntensity = new List<IsotopePair>();
        //                        for (int a = 0; a < pairs.Count(); a++)
        //                        {
        //                            otherPair = pairs[a];
        //                            for (int h = 0; h < numChannels; h += numIsotopologues)
        //                            {
        //                                for (int s = 0; s < numIsotopes; s++)
        //                                {
        //                                    otherPeak1 = otherPair.peaks[h, s];
        //                                    otherPeak2 = otherPair.peaks[h + 1, s];

        //                                    //Both channels present
        //                                    if (otherPeak1 != null && otherPeak2 != null)
        //                                    {
        //                                        isotopePair = new IsotopePair(otherPair, s, otherPeak1.X, otherPeak1.GetDenormalizedIntensity(injectionTime), otherPeak2.X, otherPeak2.GetDenormalizedIntensity(injectionTime));
        //                                        if (isotopePair.intensity <= singlePeakIntensity)
        //                                        {
        //                                            lowerIntensity.Add(isotopePair);
        //                                        }
        //                                        else
        //                                        {
        //                                            higherIntensity.Add(isotopePair);
        //                                        }
        //                                    }
        //                                    else if (otherPeak1 == null && otherPeak2 == null)
        //                                    {

        //                                    }
        //                                    //One channel present
        //                                    else
        //                                    {
        //                                        if (otherPeak1 != null)
        //                                        {
        //                                            isotopePair = new IsotopePair(otherPair, s, otherPeak1.X, otherPeak1.GetDenormalizedIntensity(injectionTime), 0, 0);
        //                                            if (isotopePair.intensity <= singlePeakIntensity)
        //                                            {
        //                                                lowerIntensity.Add(isotopePair);
        //                                            }
        //                                            else
        //                                            {
        //                                                higherIntensity.Add(isotopePair);
        //                                            }
        //                                        }
        //                                        if (otherPeak2 != null)
        //                                        {
        //                                            isotopePair = new IsotopePair(otherPair, s, otherPeak2.X, otherPeak2.GetDenormalizedIntensity(injectionTime), 0, 0);
        //                                            if (isotopePair.intensity <= singlePeakIntensity)
        //                                            {
        //                                                lowerIntensity.Add(isotopePair);
        //                                            }
        //                                            else
        //                                            {
        //                                                higherIntensity.Add(isotopePair);
        //                                            }
        //                                        }
        //                                    }
        //                                }
        //                            }
        //                        }

        //                        // Calculate the missing channel frequencies for pairs both lower and higher in intensity than the single peak
        //                        lowerFrequency = 0;
        //                        int lowerCount = 0;
        //                        for (int l = 0; l < lowerIntensity.Count(); l++)
        //                        {
        //                            if (lowerIntensity[l].missingChannel)
        //                            {
        //                                lowerCount++;
        //                            }
        //                        }
        //                        lowerFrequency = (double)lowerCount / (double)lowerIntensity.Count();

        //                        higherFrequency = 0;
        //                        int higherCount = 0;
        //                        for (int l = 0; l < higherIntensity.Count(); l++)
        //                        {
        //                            if (higherIntensity[l].missingChannel)
        //                            {
        //                                higherCount++;
        //                            }
        //                        }
        //                        higherFrequency = (double)higherCount / (double)higherIntensity.Count();

        //                        // A peak is deemed coalesced if the frequency of a missing channel is at least 1.5-fold greater for the more intense pairs than the less intense (each category must have at least 2 contributing pairs)
        //                        if (lowerIntensity.Count() > 1 && higherIntensity.Count() > 1 && higherFrequency / lowerFrequency > 1.5)
        //                        {
        //                            coalescenceDetected = true;
        //                            if (coalescedPeakIntensities == null)
        //                            {
        //                                coalescedPeakIntensities = new List<double>();
        //                            }
        //                            else
        //                            {
        //                                // Keep track of peak intensities that are deemed to be coalesced
        //                                coalescedPeakIntensities.Add(singlePeakIntensity);
        //                            }

        //                            // Set single peak to null
        //                            if (peak1 == null)
        //                            {
        //                                pair.peaks[c + 1, j] = null;
        //                            }
        //                            if (peak2 == null)
        //                            {
        //                                pair.peaks[c, j] = null;
        //                            }
        //                        }
        //                    }
        //                }

        //                if (numIsotopologues == 4)
        //                {
        //                    peak1 = pair.peaks[c, j];
        //                    peak2 = pair.peaks[c + 1, j];
        //                    peak3 = pair.peaks[c + 2, j];
        //                    peak4 = pair.peaks[c + 3, j];

        //                    // Do nothing for complete pairs, null pairs, or pairs below the intensity threshold
        //                    if (peak1 != null && peak2 != null && peak3 != null && peak4 != null)
        //                    {
                                
        //                    }
        //                    else if (peak1 == null && peak2 == null && peak3 == null && peak4 == null)
        //                    {
                                
        //                    }
        //                    //else if (log10MaxIntensity <= System.Math.Log10(Form1.MAXIMUMDNL))
        //                    //{
                                
        //                    //}
        //                    // Proceed with all other pairs
        //                    else 
        //                    {
        //                        ILabeledPeak otherPeak1;
        //                        ILabeledPeak otherPeak2;
        //                        ILabeledPeak otherPeak3;
        //                        ILabeledPeak otherPeak4;
        //                        double presentPeaksTotalIntensity = 0;
        //                        int presentPeaksCount = 0;
        //                        double presentPeaksIntensity = 0;

        //                        // Calculate the maximum intensity of all present peaks
        //                        if (peak1 != null)
        //                        {
        //                            if (peak1.GetDenormalizedIntensity(injectionTime) > presentPeaksIntensity)
        //                            {
        //                                presentPeaksTotalIntensity += peak1.GetDenormalizedIntensity(injectionTime);
        //                                presentPeaksCount++;
        //                            }
        //                        }
        //                        if (peak2 != null)
        //                        {
        //                            if (peak2.GetDenormalizedIntensity(injectionTime) > presentPeaksIntensity)
        //                            {
        //                                presentPeaksTotalIntensity += peak2.GetDenormalizedIntensity(injectionTime);
        //                                presentPeaksCount++;
        //                            }
        //                        }
        //                        if (peak3 != null)
        //                        {
        //                            if (peak3.GetDenormalizedIntensity(injectionTime) > presentPeaksIntensity)
        //                            {
        //                                presentPeaksTotalIntensity = peak3.GetDenormalizedIntensity(injectionTime);
        //                                presentPeaksCount++;
        //                            }
        //                        }
        //                        if (peak4 != null)
        //                        {
        //                            if (peak4.GetDenormalizedIntensity(injectionTime) > presentPeaksIntensity)
        //                            {
        //                                presentPeaksTotalIntensity = peak4.GetDenormalizedIntensity(injectionTime);
        //                                presentPeaksCount++;
        //                            }
        //                        }

        //                        presentPeaksIntensity = presentPeaksTotalIntensity / ((double)presentPeaksCount);

        //                        // Sort the peptide's other pairs (either both or one channel present) based on their intensities relative to the single peak
        //                        lowerIntensity = new List<IsotopePair>();
        //                        higherIntensity = new List<IsotopePair>();
        //                        Pair otherPair;
        //                        for (int a = 0; a < pairs.Count(); a++)
        //                        {
        //                            otherPair = pairs[a];
        //                            for (int h = 0; h < numChannels; h += numIsotopologues)
        //                            {
        //                                for (int s = 0; s < numIsotopes; s++)
        //                                {
        //                                    otherPeak1 = otherPair.peaks[h, s];
        //                                    otherPeak2 = otherPair.peaks[h + 1, s];
        //                                    otherPeak3 = otherPair.peaks[h + 2, s];
        //                                    otherPeak4 = otherPair.peaks[h + 3, s];

        //                                    // All channels present
        //                                    if (otherPeak1 != null && otherPeak2 != null && otherPeak3 != null && otherPeak4 != null)
        //                                    {
        //                                        isotopePair = new IsotopePair(otherPair, s, otherPeak1.X, otherPeak1.GetDenormalizedIntensity(injectionTime), otherPeak2.X, otherPeak2.GetDenormalizedIntensity(injectionTime), otherPeak3.X, otherPeak3.GetDenormalizedIntensity(injectionTime), otherPeak4.X, otherPeak4.GetDenormalizedIntensity(injectionTime));
        //                                        if (isotopePair.intensity <= presentPeaksIntensity)
        //                                        {
        //                                            lowerIntensity.Add(isotopePair);
        //                                        }
        //                                        else
        //                                        {
        //                                            higherIntensity.Add(isotopePair);
        //                                        }
        //                                    }
        //                                    // No channels present
        //                                    else if (otherPeak1 == null && otherPeak2 == null && otherPeak3 == null && otherPeak4 == null)
        //                                    {
        //                                        // Do nothing
        //                                    }
        //                                    // Incomplete pairs -- set null peaks to blank peaks
        //                                    else
        //                                    {
        //                                        MZPeak blank = new MZPeak(0,0);

        //                                        if (otherPeak1 == null)
        //                                        {
        //                                            otherPeak1 = (ILabeledPeak)blank;
        //                                        }
        //                                        if (otherPeak2 == null)
        //                                        {
        //                                            otherPeak2 = (ILabeledPeak)blank;
        //                                        }
        //                                        if (otherPeak3 == null)
        //                                        {
        //                                            otherPeak3 = (ILabeledPeak)blank;
        //                                        }
        //                                        if (otherPeak4 == null)
        //                                        {
        //                                            otherPeak4 = (ILabeledPeak)blank;
        //                                        }

        //                                        isotopePair = new IsotopePair(otherPair, s, otherPeak1.X, otherPeak1.GetDenormalizedIntensity(injectionTime), otherPeak2.X, otherPeak2.GetDenormalizedIntensity(injectionTime), otherPeak3.X, otherPeak3.GetDenormalizedIntensity(injectionTime), otherPeak4.X, otherPeak4.GetDenormalizedIntensity(injectionTime));
        //                                        if (isotopePair.intensity <= presentPeaksIntensity)
        //                                        {
        //                                            lowerIntensity.Add(isotopePair);
        //                                        }
        //                                        else
        //                                        {
        //                                            higherIntensity.Add(isotopePair);
        //                                        }
        //                                    }
        //                                }
        //                            }
        //                        }

        //                        // Calculate the missing channel frequencies for pairs both lower and higher in intensity than the single peak
        //                        lowerFrequency = 0;
        //                        int lowerCount = 0;
        //                        for (int l = 0; l < lowerIntensity.Count(); l++)
        //                        {
        //                            if (lowerIntensity[l].missingChannel)
        //                            {
        //                                lowerCount++;
        //                            }
        //                        }
        //                        lowerFrequency = (double)lowerCount / (double)lowerIntensity.Count();

        //                        higherFrequency = 0;
        //                        int higherCount = 0;
        //                        for (int l = 0; l < higherIntensity.Count(); l++)
        //                        {
        //                            if (higherIntensity[l].missingChannel)
        //                            {
        //                                higherCount++;
        //                            }
        //                        }
        //                        higherFrequency = (double)higherCount / (double)higherIntensity.Count();

        //                        // A pair is deemed coalesced if the frequency of a missing channel is at least 1.5-fold greater for the more intense pairs than the less intense (each category must have at least 2 contributing pairs)
        //                        if (lowerIntensity.Count() > 1 && higherIntensity.Count() > 1 && higherFrequency / lowerFrequency > 1.5)
        //                        {
        //                            coalescenceDetected = true;
        //                            if (coalescedPeakIntensities == null)
        //                            {
        //                                coalescedPeakIntensities = new List<double>();
        //                            }
        //                            else
        //                            {
        //                                coalescedPeakIntensities.Add(presentPeaksIntensity);
        //                            }

        //                            // Set single peaks to null for a coalesced pair
        //                            pair.peaks[c, j] = null;
        //                            pair.peaks[c + 1, j] = null;
        //                            pair.peaks[c + 2, j] = null;
        //                            pair.peaks[c + 3, j] = null;
        //                        }
        //                    }
        //                }
        //            } // End channel loop
        //        } // End isotope loop
        //    } // End pair loop
        //}

        /* Checks the pair spacing to make sure it falls within the spacing tolerance
         */
        public void checkPairSpacing(MSDataFile rawFile, List<Spacing> spacings = null)
        {            
            int charge = bestPSMs[rawFile].Charge;
            List<Pair> pairs;
            allPairs.TryGetValue(rawFile, out pairs);
            int peakCount;
            bool[,] spacingChecks;
            bool spacingChecked;

            if (numIsotopologues < 2 && clusterLabel)
            {
                foreach (Pair pair in pairs)
                {
                    foreach (IsotopePair isotopePair in pair.IsotopePairs)
                    {
                        peakCount = isotopePair.TotalPeakCount;
                        spacingChecks = isotopePair.GetSpacingCheckArray();

                        for (int j = 0; j < numChannels; j++)
                        {
                            spacingChecked = false;
                            for (int c = 0; c < numChannels; c++)
                            {
                                if (spacingChecks[j, c]) spacingChecked = true;
                            }
                            if (isotopePair.ChannelPeaks[j] != null && !spacingChecked)
                            {
                                peakCount--;
                                isotopePair.ChannelPeaks[j] = null;
                            }
                        }

                        if (peakCount < peaksNeeded)
                        {
                            for (int j = 0; j < numChannels; j++)
                            {
                                isotopePair.ChannelPeaks[j] = null;
                            }
                        }
                    }
                }
            }
            else if (numIsotopologues > 1 && isotopologueLabel)
            {
                foreach (Pair pair in pairs)
                {
                    foreach (IsotopePair isotopePair in pair.IsotopePairs)
                    {
                        for (int c = 0; c < numClusters; c++)
                        {
                            peakCount = isotopePair.PeakCountByCluster[c];
                            spacingChecks = isotopePair.GetSpacingCheckArray(c);
                            int channelIndex = c * numIsotopologues;

                            for (int j = 0; j < numIsotopologues; j++)
                            {
                                spacingChecked = false;
                                for (int m = 0; m < numIsotopologues; m++)
                                {
                                    if (spacingChecks[j, m]) spacingChecked = true;
                                }
                                if (isotopePair.ChannelPeaks[j + channelIndex] != null && !spacingChecked)
                                {
                                    peakCount--;
                                    isotopePair.ChannelPeaks[j + channelIndex] = null;
                                }
                            }

                            if (peakCount < peaksNeeded)
                            {
                                for (int j = channelIndex; j < channelIndex + numIsotopologues; j++)
                                {
                                    isotopePair.ChannelPeaks[j] = null;
                                }
                            }
                        }
                    }
                }
            }

            List<Pair> spacingsCheckedList = new List<Pair>();
            List<Pair> completePairList = null;
            completePairs.TryGetValue(rawFile, out completePairList);

            if (numIsotopologues > 1 && numClusters > 1)
            {
                foreach (Pair currentPair in pairs)
                {
                    foreach (IsotopePair currentIsotopePair in currentPair.IsotopePairs)
                    {
                        for (int c = 0; c < numClusters; c++)
                        {
                            if (currentIsotopePair.CompleteByCluster[c])
                            {
                                completePairList = addPairToList(completePairList, currentIsotopePair, c, rawFile);
                            }
                            else if (currentIsotopePair.PeakCountByCluster[c] >= peaksNeeded)
                            {
                                spacingsCheckedList = addPairToList(spacingsCheckedList, currentIsotopePair, c, rawFile);
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
                    foreach (IsotopePair currentIsotopePair in currentPair.IsotopePairs)
                    {
                        if (currentIsotopePair.Complete)
                        {
                            completePairList = addPairToList(completePairList, currentIsotopePair, 0, rawFile);
                        }
                        else if (currentIsotopePair.TotalPeakCount >= peaksNeeded)
                        {
                            spacingsCheckedList = addPairToList(spacingsCheckedList, currentIsotopePair, 0, rawFile);
                        }
                    }
                }
                completePairs[rawFile] = completePairList;
                allPairs[rawFile] = spacingsCheckedList;
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

        //public void checkIsotopeDistribution(MSDataFile rawFile, List<Pair> pairs)
        //{
        //    List<Pair> goodIsotopeDistributionPairs = new List<Pair>();

        //    if (numIsotopologues > 1)
        //    {
        //        foreach (Pair pair in pairs)
        //        {
        //            for (int c = 0; c < numClusters; c++)
        //            {
        //                if (numIsotopes > 1 && pair.checkIsotopeDistribution(c))
        //                {
        //                    for (int j = 0; j < numIsotopes; j++)
        //                    {
        //                        if (pair.peakCount[c, j] < peaksNeeded)
        //                        {
        //                            int channelIndex = c * numIsotopologues;
        //                            for (int i = channelIndex; i < channelIndex + numIsotopologues; i++)
        //                            {
        //                                pair.peaks[i, j] = null;
        //                            }
        //                        }
        //                    }
        //                }
        //            }
        //            if (pair.totalPeakCount > 0)
        //            {
        //                goodIsotopeDistributionPairs.Add(pair);
        //            }
        //        }
        //    }
        //    else
        //    {
        //        foreach (Pair pair in pairs)
        //        {
        //            for (int c = 0; c < numChannels; c++)
        //            {
        //                pair.checkIsotopeDistribution(c);
        //            }

        //            for (int j = 0; j < numIsotopes; j++)
        //            {
        //                if (pair.peakCount[0, j] < peaksNeeded)
        //                {
        //                    for (int c = 0; c < numChannels; c++)
        //                    {
        //                        pair.peaks[c, j] = null;
        //                    }
        //                }
        //            }

        //            if (pair.totalPeakCount > 0)
        //            {
        //                goodIsotopeDistributionPairs.Add(pair);
        //            }
        //        }
        //    }
        //    allPairs[rawFile] = goodIsotopeDistributionPairs;
        //}

        public void sortPairs(MSDataFile rawFile)
        {
            PeakPattern[,] pairPatterns;
            List<Pair> associatedPairs = null;

            if (allPairs.TryGetValue(rawFile, out associatedPairs))
            {
                // NeuCode quantification
                if (numIsotopologues > 1 && numIsotopologues < 6)
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
                                foreach (IsotopePair isotopePair in pair.IsotopePairs)
                                {
                                    if (isotopePair.PeakCountByCluster[c] >= peaksNeeded && isotopePair.PeakCountByCluster[c] < numIsotopologues && isotopePair.PeakPattern[c] >= 0)
                                    {
                                        IsotopePair cleanedPair = new IsotopePair(isotopePair.Parent, isotopePair.Isotope);
                                        for (int i = 0; i < numChannels; i++)
                                        {
                                            if (i >= channelIndex && i < channelIndex + numIsotopologues)
                                            {
                                                cleanedPair.ChannelPeaks[i] = isotopePair.ChannelPeaks[i];
                                            }
                                        }
                                        pairPatterns[isotopePair.PeakCountByCluster[c], isotopePair.PeakPattern[c]].PatternPairs.Add(cleanedPair);
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

                                foreach (IsotopePair isotopePair in pairPatterns[bestNumPeaks, bestPatternPossibilities].PatternPairs)
                                {
                                    if (!newAllPairs.Contains(isotopePair.Parent))
                                    {
                                        Pair pairToAdd = new Pair(this, rawFile, isotopePair.Parent.ScanNumber, isotopePair.Parent.Charge, isotopePair.Parent.InjectionTime, isotopePair.Parent.RetentionTime);
                                        pairToAdd.IsotopePairs.Add(isotopePair);
                                        newAllPairs.Add(pairToAdd);
                                    }
                                    else
                                    {
                                        newAllPairs.Find(pair => pair.ScanNumber == isotopePair.Parent.ScanNumber).IsotopePairs.Add(isotopePair);
                                    }
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
                else if (numIsotopologues < 2)
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
                            foreach (IsotopePair isotopePair in pair.IsotopePairs)
                            {
                                if (isotopePair.PeakCountByCluster[0] >= peaksNeeded && isotopePair.PeakCountByCluster[0] < numIsotopologues && isotopePair.PeakPattern[0] >= 0)
                                {
                                    IsotopePair cleanedPair = new IsotopePair(isotopePair.Parent, isotopePair.Isotope);
                                    for (int i = 0; i < numChannels; i++)
                                    {
                                        cleanedPair.ChannelPeaks[i] = isotopePair.ChannelPeaks[i];
                                    }
                                    pairPatterns[isotopePair.PeakCountByCluster[0], isotopePair.PeakPattern[0]].PatternPairs.Add(cleanedPair);
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
                        double secondBestScore = 0;
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
                            missingChannelPattern[0, 0] = bestNumPeaks;
                            missingChannelPattern[0, 1] = bestPatternPossibilities;

                            foreach (IsotopePair isotopePair in pairPatterns[bestNumPeaks, bestPatternPossibilities].PatternPairs)
                            {
                                if (!newAllPairs.Contains(isotopePair.Parent))
                                {
                                    Pair pairToAdd = new Pair(this, rawFile, isotopePair.Parent.ScanNumber, isotopePair.Parent.Charge, isotopePair.Parent.InjectionTime, isotopePair.Parent.RetentionTime);
                                    pairToAdd.IsotopePairs.Add(isotopePair);
                                    newAllPairs.Add(pairToAdd);
                                }
                                else
                                {
                                    newAllPairs.Find(pair => pair.ScanNumber == isotopePair.Parent.ScanNumber).IsotopePairs.Add(isotopePair);
                                }
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

        public void setTheoreticalIsotopicDistribution()
        {
            theoreticalIsotopicDistribution = new double[numIsotopes];

            int countC12 = Pair.countElement("C", sequenceNoMods);
            int countH1 = Pair.countElement("H", sequenceNoMods);
            int countN14 = Pair.countElement("N", sequenceNoMods);
            int countO16 = Pair.countElement("O", sequenceNoMods);
            int countS32 = Pair.countElement("S", sequenceNoMods);
            double[] relativeAbundances = new double[numIsotopes];
    
            if (numIsotopes == 3)
            {  
                relativeAbundances[0] = 100.0;
                relativeAbundances[1] = (countC12 * 1.1) + (countH1 * 0.015) + (countN14 * 0.37);
                relativeAbundances[2] = (Math.Pow(countC12 * 1.1, 2)) / 200.0 + (countO16 * 0.2) + (countS32 * 4.21);
            }

            for (int i = 0; i < numIsotopes; i++)
            {
                theoreticalIsotopicDistribution[i] = relativeAbundances[i];
            }
        }

        public List<Pair> addPairToList(List<Pair> initialList, IsotopePair isotopePair, int cluster, MSDataFile rawFile)
        {
            List<Pair> updatedList = new List<Pair>();
            if (initialList == null) initialList = new List<Pair>();
            updatedList = initialList;
            Pair pairToAdd = new Pair(this, rawFile, isotopePair.Parent.ScanNumber, isotopePair.Parent.Charge, isotopePair.Parent.InjectionTime, isotopePair.Parent.RetentionTime);

            if (numIsotopologues > 1 && numClusters > 1)
            {
                if (updatedList.Count == 0)
                {
                    pairToAdd.IsotopePairs.Add(isotopePair);
                    updatedList.Add(pairToAdd);
                }
                else
                {
                    bool found = false;
                    foreach (Pair oldPair in updatedList)
                    {
                        if (oldPair.ScanNumber == pairToAdd.ScanNumber)
                        {
                            bool isotopeFound = false;
                            foreach (IsotopePair oldIsotopePair in oldPair.IsotopePairs)
                            {
                                if (oldIsotopePair.Isotope == isotopePair.Isotope)
                                {
                                    isotopeFound = true;
                                    int channelIndex = cluster * numIsotopologues;
                                    for (int ch = channelIndex; ch < channelIndex + numIsotopologues; ch++)
                                    {
                                        oldIsotopePair.ChannelPeaks[ch] = isotopePair.ChannelPeaks[ch];
                                    }
                                }
                            }
                            if (!isotopeFound) oldPair.IsotopePairs.Add(isotopePair);
                            found = true;
                        }
                    }
                    if (!found)
                    {
                        pairToAdd.IsotopePairs.Add(isotopePair);
                        updatedList.Add(pairToAdd);
                    }
                }
            }
            else
            {
                if (updatedList.Count == 0)
                {
                    pairToAdd.IsotopePairs.Add(isotopePair);
                    updatedList.Add(pairToAdd);
                }
                else
                {
                    bool found = false;
                    foreach (Pair oldPair in updatedList)
                    {
                        if (oldPair.ScanNumber == pairToAdd.ScanNumber)
                        {
                            oldPair.IsotopePairs.Add(isotopePair);
                            found = true;
                        }
                    }
                    if (!found)
                    {
                        pairToAdd.IsotopePairs.Add(isotopePair);
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
            
            //Apply noise to missing channels
            List<Pair> all;
            List<Pair> completeOnly;
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
                    foreach (IsotopePair isotopePair in pair.IsotopePairs)
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
                            if (!isotopePair.CompleteByCluster[c] && isotopePair.PeakCountByCluster[c] >= peaksNeeded)
                            {
                                // Check for coalescence
                                if (coalescenceDetected && isotopePair.MaxIntensityByCluster[c] > coalescenceIntensity)
                                {
                                    for (int m = channelIndex; m < channelIndex + numIsotopologues; m++)
                                    {
                                        isotopePair.ChannelPeaks[m] = null;
                                    }
                                }
                                // Apply noise to missing channels
                                else
                                {
                                    for (int m = channelIndex; m < channelIndex + numIsotopologues; m++)
                                    {
                                        if (isotopePair.ChannelPeaks[m] == null)
                                        {
                                            ThermoLabeledPeak noisePeak = new ThermoLabeledPeak(0.0, pair.AverageNoise, pair.Charge, pair.AverageNoise);
                                            isotopePair.ChannelPeaks[m] = noisePeak;
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
                    foreach (IsotopePair isotopePair in pair.IsotopePairs)
                    {
                        // For incomplete pairs
                        if (!isotopePair.Complete && isotopePair.TotalPeakCount >= peaksNeeded)
                        {
                            for (int m = 0; m < numChannels; m++)
                            {
                                if (isotopePair.ChannelPeaks[m] == null)
                                {
                                    ThermoLabeledPeak noisePeak = new ThermoLabeledPeak(0.0, pair.AverageNoise, pair.Charge, pair.AverageNoise);
                                    isotopePair.ChannelPeaks[m] = noisePeak;
                                }
                            }
                        }
                        
                        //// For complete pairs
                        //if (pair.complete[0, j])
                        //{
                        //    // First check to see if the peptide is already on the complete list
                        //    bool found = false;
                        //    int SN = pair.scanNumber;

                        //    foreach (Pair noNBC in completeOnly)
                        //    {
                        //        if (noNBC.scanNumber == SN)
                        //        {
                        //            // If found, update peaks for that pair
                        //            for (int m = 0; m < numChannels; m++)
                        //            {
                        //                noNBC.peaks[m, j] = pair.peaks[m, j];
                        //            }
                        //            found = true;
                        //        }
                        //    }

                        //    // If not found, add pair
                        //    if (!found)
                        //    {
                        //        Pair noNBCPair = new Pair(this, rawFile, pair.scanNumber, pair.injectionTime, pair.retentionTime);
                        //        for (int m = 0; m < numChannels; m++)
                        //        {
                        //            noNBCPair.peaks[m, j] = pair.peaks[m, j];
                        //        }
                        //        completeOnly.Add(noNBCPair);
                        //    }
                        //}

                        //// For incomplete pairs
                        //else if (pair.peakCount[0, j] >= peaksNeeded)
                        //{
                        //    // Apply noise to missing channels
                        //    for (int m = 0; m < numChannels; m++)
                        //    {
                        //        if (pair.peaks[m, j] == null)
                        //        {
                        //            ThermoLabeledPeak noisePeak = new ThermoLabeledPeak(0.0, pair.averageNoise, pair.charge, pair.averageNoise);
                        //            pair.peaks[m, j] = noisePeak;
                        //        }
                        //    }
                        //}
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
            if (numIsotopologues > 1 && isotopologueLabel)
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
                int numPeptideVersions = numChannels;
                if (!clusterLabel) numPeptideVersions = 1;
                ILabeledPeak[] peaks = new ILabeledPeak[numPeptideVersions];
                int peaksCount = 0;
                
                //if (Form1.SILAC_DUPLEX_LEUCN && Form1.CORRECTLEUNLOSS)
                //{
                //    int correction = correctIsotopeDistributions(precursorScan, rawFile);
                //    theoMasses[1, 0] -= correction * (Constants.NITROGEN15 - Constants.NITROGEN);
                //}
                for (int i = 0; i < numPeptideVersions; i++)
                {
                    peaks[i] = largestPeak(theoMasses[i, 0], precursorScan, firstSearchMassRange, rawFile);
                    if (peaks[i] != null) peaksCount++;
                }
                if (peaksCount == numPeptideVersions)
                {
                    for (int i = 0; i < numPeptideVersions; i++)
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
            if (!labeled) return false;
            else if (numIsotopologues < 2 && clusterLabel) return true;
            else if (numIsotopologues > 1 && !isotopologueLabel) return false;
            else
            {
                int charge = bestPSM.Charge;
                double experimentalSeparation = (spacingMassRange[1, 0].Mean) / (double)charge;
                int clusterIndex = cluster * numIsotopologues;
                double coefficient = (Math.Sqrt(2 * Math.Log(100.0 / quantSeparation))) / (Math.Sqrt(2 * Math.Log(2)));
                double theoreticalSeparation = coefficient * ((Mass.MzFromMass(theoMasses[clusterIndex, 0], charge) / (quantResolution * Math.Sqrt(400 / Mass.MzFromMass(theoMasses[clusterIndex, 0], charge)))));

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
                    foreach (IsotopePair isotopePair in pair.IsotopePairs)
                    {
                        for (int c = 0; c < numClusters; c++)
                        {
                            // Complete pairs
                            if (isotopePair.CompleteByCluster[c])
                            {
                                pairTotalCount++;
                            }
                            // Pairs with fewer than half the total number of isotopologues
                            else if (Form1.NOISEBANDCAP && isotopePair.PeakCountByCluster[c] < peaksNeeded)
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
        public bool quantFilter(IsotopePair isotopePair, int cluster, bool complete)
        {
            // Find the maximum intensity for each channel, considering incomplete sets only in the absence of complete pairs
            double[,] max = absoluteMaxIntensity;
            
            // Set maximum as 0 for isotope quant
            if (Form1.ISOTOPEQUANT)
            {
                max = new double[numChannels, 2];
            }
            
            int channelIndex;

            // Should not have any null peaks in the pair at this point
            if (numIsotopologues > 1 && isotopePair.PeakCountByCluster[cluster] != numIsotopologues)
            {
                return false;
            }
            else if (numIsotopologues < 2 && isotopePair.TotalPeakCount != numChannels)
            {
                return false;
            }

            // If there are more than 2 complete pairs, use their max intensity
            //if (countCompleteIsotopes[cluster] > peaksNeeded)
            //{
            //    max = maxCompleteIntensity;
            //}
            //else
            //{
            //    max = maxIntensity;
            //}

            //if (maxCompleteIntensity[0, 0] > 0.0 || maxCompleteIntensity[1,0] > 0.0) max = maxCompleteIntensity;
            //else max = maxIntensity;

            bool includePair = false;

            // NeuCode + multiple clusters
            if (numIsotopologues > 1 && numClusters > 1)
            {
                channelIndex = cluster * numIsotopologues;

                // For complete pairs
                if (isotopePair.CompleteByCluster[cluster])
                {
                    for (int i = channelIndex; i < channelIndex + numIsotopologues; i++)
                    {
                        double intensityThreshold = INTENSITYCUTOFF * max[i, 0];
                        double peakIntensity = isotopePair.ChannelPeaks[i].GetDenormalizedIntensity(isotopePair.Parent.InjectionTime);
                        if (peakIntensity >= intensityThreshold)
                        {
                            includePair = true;
                            if (complete)
                            {
                                completeTotalIntensity[i, isotopePair.Isotope] += peakIntensity;
                                completeTotalIntensity[i, numIsotopes] += peakIntensity;
                            }
                            else
                            {
                                totalIntensity[i, isotopePair.Isotope] += peakIntensity;
                                totalIntensity[i, numIsotopes] += peakIntensity;
                            }
                        }
                    }
                }
                
                // For incomplete pairs, only consider non noise-band capped channels for quantitative filtering
                else
                {
                    for (int i = channelIndex; i < channelIndex + numIsotopologues; i++)
                    {
                        double intensityThreshold = INTENSITYCUTOFF * max[i, 0];
                        double peakIntensity = isotopePair.ChannelPeaks[i].GetDenormalizedIntensity(isotopePair.Parent.InjectionTime);                        
                        if (peakIntensity >= intensityThreshold)
                        {
                            includePair = true;
                            if (complete)
                            {
                                completeTotalIntensity[i, isotopePair.Isotope] += peakIntensity;
                                completeTotalIntensity[i, numIsotopes] += peakIntensity;
                            }
                            else
                            {
                                totalIntensity[i, isotopePair.Isotope] += peakIntensity;
                                totalIntensity[i, numIsotopes] += peakIntensity;
                            }
                        }
                    }
                }
            }
            // Single-cluster NeuCode and traditional SILAC
            else
            {
                for (int i = 0; i < numChannels; i++)
                {
                    // For complete pairs
                    if (isotopePair.Complete)
                    {
                        double intensityThreshold = INTENSITYCUTOFF * max[i, 0];
                        double peakIntensity = isotopePair.ChannelPeaks[i].GetDenormalizedIntensity(isotopePair.Parent.InjectionTime);
                        if (peakIntensity >= intensityThreshold)
                        {
                            includePair = true;
                            if (complete)
                            {
                                completeTotalIntensity[i, isotopePair.Isotope] += peakIntensity;
                                completeTotalIntensity[i, numIsotopes] += peakIntensity;
                            }
                            else
                            {
                                totalIntensity[i, isotopePair.Isotope] += peakIntensity;
                                totalIntensity[i, numIsotopes] += peakIntensity;
                            }
                        }
                    }
                    // For incomplete pairs, only consider non noise-band capped channels for quantitative filtering
                    else
                    {
                        double intensityThreshold = INTENSITYCUTOFF * max[i, 0];
                        double peakIntensity = isotopePair.ChannelPeaks[i].GetDenormalizedIntensity(isotopePair.Parent.InjectionTime);
                        if (peakIntensity >= intensityThreshold)
                        {
                            includePair = true;
                            if (complete)
                            {
                                completeTotalIntensity[i, isotopePair.Isotope] += peakIntensity;
                                completeTotalIntensity[i, numIsotopes] += peakIntensity;
                            }
                            else
                            {
                                totalIntensity[i, isotopePair.Isotope] += peakIntensity;
                                totalIntensity[i, numIsotopes] += peakIntensity;
                            }
                        }
                    }                    
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
            return includePair;
        }

        public void filterPairs(MSDataFile rawFile)
        {
            List<Pair> unfilteredPairs = allPairs[rawFile];
            List<Pair> rTFilteredPairs = new List<Pair>();
            List<Pair> monoisotopeFilteredPairs = new List<Pair>();
            List<Pair> isotopicDistributionFilteredPairs = new List<Pair>();
            List<MSDataScan> MS1Scans = fullScanList;
            List<int> MS1ScanNumbers = new List<int>();
            foreach (MSDataScan Scan in MS1Scans)
            {
                MS1ScanNumbers.Add(Scan.SpectrumNumber);
            }
            MS1ScanNumbers.Sort();
            List<int> PSMScanNumbers = new List<int>();
            foreach (PeptideSpectralMatch Psm in PSMs[rawFile])
            {
                PSMScanNumbers.Add(Psm.MS1ScanNumber);
            }
            PSMScanNumbers.Sort();
            List<int> MS1ScanNumbersBeforePSM = new List<int>();
            List<Pair> PairsBeforePSM = new List<Pair>();
            List<int> MS1ScanNumbersAfterPSM = new List<int>();
            List<Pair> PairsAfterPSM = new List<Pair>();
            List<int> MS1ScanNumbersContainingPSM = new List<int>();
            List<Pair> PairsContainingPSM = new List<Pair>();
            int firstPSMScanNumber = PSMScanNumbers[0];
            int lastPSMScanNumber = PSMScanNumbers[PSMScanNumbers.Count - 1];
            if (lastPSMScanNumber - firstPSMScanNumber >= 1000)
            {
                PSMScanNumbers.Clear();
                PSMScanNumbers.Add(bestPSMs[rawFile].MS1ScanNumber);
                firstPSMScanNumber = bestPSMs[rawFile].MS1ScanNumber;
                lastPSMScanNumber = firstPSMScanNumber;
            }
            foreach (int Scan in MS1ScanNumbers)
            {
                if (Scan < firstPSMScanNumber) MS1ScanNumbersBeforePSM.Add(Scan);
                else if (Scan > lastPSMScanNumber) MS1ScanNumbersAfterPSM.Add(Scan);
                else MS1ScanNumbersContainingPSM.Add(Scan);
            }
            foreach (Pair pair in unfilteredPairs)
            {
                if (pair.ScanNumber < firstPSMScanNumber) PairsBeforePSM.Add(pair);
                else if (pair.ScanNumber > lastPSMScanNumber) PairsAfterPSM.Add(pair);
                else PairsContainingPSM.Add(pair);
            }

            // Sort scan & pairs
            List<int> SortedMS1ScanNumbersBeforePSM = new List<int>();
            for (int i = MS1ScanNumbersBeforePSM.Count - 1; i >= 0; i--)
            {
                SortedMS1ScanNumbersBeforePSM.Add(MS1ScanNumbersBeforePSM[i]);
            }
            List<Pair> SortedPairsBeforePSM = PairsBeforePSM.OrderByDescending(pair => pair.ScanNumber).ToList();

            foreach (Pair pair in PairsContainingPSM)
            {
                rTFilteredPairs.Add(pair);
            }


            int currentScanIndex;
            int previousScanIndex;
            Pair currentPair;

            bool stopIncludingPairs = false;
            int currentIndex = 0;
            int numSkippedScans = 0;
            while (!stopIncludingPairs && currentIndex < SortedPairsBeforePSM.Count)
            {
                currentPair = SortedPairsBeforePSM[currentIndex];
                currentScanIndex = SortedMS1ScanNumbersBeforePSM.IndexOf(currentPair.ScanNumber);
                
                if (currentIndex == 0)
                {
                    numSkippedScans = currentScanIndex;
                }
                else
                {
                    previousScanIndex = SortedMS1ScanNumbersBeforePSM.IndexOf(SortedPairsBeforePSM[currentIndex - 1].ScanNumber);
                    numSkippedScans = currentScanIndex - previousScanIndex;
                }

                if (numSkippedScans > 1) stopIncludingPairs = true;
                else
                {
                    rTFilteredPairs.Add(currentPair);
                }
                currentIndex++;
            }

            stopIncludingPairs = false;
            currentIndex = 0;
            numSkippedScans = 0;
            while (!stopIncludingPairs && currentIndex < PairsAfterPSM.Count)
            {
                currentPair = PairsAfterPSM[currentIndex];
                currentScanIndex = MS1ScanNumbersAfterPSM.IndexOf(currentPair.ScanNumber);

                if (currentIndex == 0)
                {
                    numSkippedScans = currentScanIndex;
                }
                else
                {
                    previousScanIndex = MS1ScanNumbersAfterPSM.IndexOf(PairsAfterPSM[currentIndex - 1].ScanNumber);
                    numSkippedScans = currentScanIndex - previousScanIndex;
                }
                if (numSkippedScans > 1) stopIncludingPairs = true;
                else
                {
                    rTFilteredPairs.Add(currentPair);
                }
                currentIndex++;
            }

            //for (int i = 0; i < unfilteredPairs.Count; i++)
            //{
            //    Pair currentPair = unfilteredPairs[i];
            //    Pair previousPair;
            //    Pair nextPair;
            //    int currentPairIndex = MS1ScanNumbers.IndexOf(currentPair.ScanNumber);
            //    int nextPairIndex;
            //    int previousPairIndex;


            //    if (unfilteredPairs.Count == 1) rTFilteredPairs.Add(currentPair);
            //    else if (i == 0)
            //    {
            //        nextPair = unfilteredPairs[i + 1];
            //        nextPairIndex = MS1ScanNumbers.IndexOf(nextPair.ScanNumber);

            //        if (nextPairIndex - currentPairIndex <= 2 || currentPair.PeaksPerIsotope[0] > 0) rTFilteredPairs.Add(currentPair);

            //    }
            //    else if (i == unfilteredPairs.Count - 1)
            //    {
            //        previousPair = unfilteredPairs[i - 1];
            //        previousPairIndex = MS1ScanNumbers.IndexOf(previousPair.ScanNumber);

            //        if (currentPairIndex - previousPairIndex <= 2 || currentPair.PeaksPerIsotope[0] > 0) rTFilteredPairs.Add(currentPair);
            //    }
            //    else
            //    {
            //        nextPair = unfilteredPairs[i + 1];
            //        nextPairIndex = MS1ScanNumbers.IndexOf(nextPair.ScanNumber);
            //        previousPair = unfilteredPairs[i - 1];
            //        previousPairIndex = MS1ScanNumbers.IndexOf(previousPair.ScanNumber);

            //        if (nextPairIndex - currentPairIndex <= 2 || currentPairIndex - previousPairIndex <= 2 || currentPair.PeaksPerIsotope[0] > 0) rTFilteredPairs.Add(currentPair);
            //    }
            //}

            // Sort retention time filtered pairs
            List<Pair> sortedRTFilteredPairs = rTFilteredPairs.OrderBy(pair => pair.ScanNumber).ToList();
            
            if (sortedRTFilteredPairs.Count <= 1)
            {
                monoisotopeFilteredPairs = sortedRTFilteredPairs;
            }
            else
            {
                int firstMonoScanNumber = sortedRTFilteredPairs[sortedRTFilteredPairs.Count - 1].ScanNumber;
                int lastMonoScanNumber = sortedRTFilteredPairs[0].ScanNumber;

                for (int i = 0; i < sortedRTFilteredPairs.Count; i++)
                {
                    if (sortedRTFilteredPairs[i].PeaksPerIsotope[0] > 0 && sortedRTFilteredPairs[i].ScanNumber < firstMonoScanNumber) firstMonoScanNumber = sortedRTFilteredPairs[i].ScanNumber;
                    if (sortedRTFilteredPairs[i].PeaksPerIsotope[0] > 0 && sortedRTFilteredPairs[i].ScanNumber > lastMonoScanNumber) lastMonoScanNumber = sortedRTFilteredPairs[i].ScanNumber;
                }

                if (firstMonoScanNumber <= lastMonoScanNumber)
                {
                    Range<int> scanNumberRange = new Range<int>(firstMonoScanNumber, lastMonoScanNumber);
                    for (int i = 0; i < sortedRTFilteredPairs.Count; i++)
                    {
                        if (scanNumberRange.Contains(sortedRTFilteredPairs[i].ScanNumber)) monoisotopeFilteredPairs.Add(sortedRTFilteredPairs[i]);
                    }
                }
                else
                {
                    monoisotopeFilteredPairs = sortedRTFilteredPairs;
                }
            }

            if (monoisotopeFilteredPairs.Count <= 1 || numIsotopes < 2)
            {
                isotopicDistributionFilteredPairs = monoisotopeFilteredPairs;
            }
            else
            {
                foreach (Pair pair in monoisotopeFilteredPairs)
                {
                    if (pair.PeaksPerIsotope[0] > 0 && pair.IsotopePairs.Count > 1)
                    {
                        List<IsotopePair> cleanedIsotopePairs = new List<IsotopePair>();
                        double monoMaxIntensity = pair.IsotopePairs[0].MaxIntensity;
                        if (pair.IsotopePairs[0].Complete)
                        {           
                            cleanedIsotopePairs.Add(pair.IsotopePairs[0]);
                            for (int i = 1; i < pair.IsotopePairs.Count; i++)
                            {
                                if (pair.IsotopePairs[i].Complete)
                                {
                                    double distributionRatio = theoreticalIsotopicDistribution[i] / theoreticalIsotopicDistribution[0];
                                    if (distributionRatio < 1)
                                    {
                                        if ((pair.IsotopePairs[i].MaxIntensity / monoMaxIntensity) <= distributionRatio * 2.0)
                                        {
                                            cleanedIsotopePairs.Add(pair.IsotopePairs[i]);
                                        }
                                    }
                                    else
                                    {
                                        if ((pair.IsotopePairs[i].MaxIntensity / monoMaxIntensity) >= distributionRatio * 0.5)
                                        {
                                            cleanedIsotopePairs.Add(pair.IsotopePairs[i]);
                                        }
                                    }
                                }
                            }
                        }
                        else
                        {
                            bool completeIsotopeFound = false;
                            for (int i = 1; i < pair.IsotopePairs.Count; i++)
                            {
                                double distributionRatio = theoreticalIsotopicDistribution[i] / theoreticalIsotopicDistribution[0];
                                if (distributionRatio < 1)
                                {
                                    if ((pair.IsotopePairs[i].MaxIntensity / monoMaxIntensity) <= distributionRatio * 2.0)
                                    {
                                        if (pair.IsotopePairs[i].Complete) completeIsotopeFound = true;
                                        cleanedIsotopePairs.Add(pair.IsotopePairs[i]);
                                    }
                                }
                                else
                                {
                                    if ((pair.IsotopePairs[i].MaxIntensity / monoMaxIntensity) >= distributionRatio * 0.5)
                                    {
                                        if (pair.IsotopePairs[i].Complete) completeIsotopeFound = true;
                                        cleanedIsotopePairs.Add(pair.IsotopePairs[i]);
                                    }
                                }
                            }

                            if (!completeIsotopeFound)
                            {
                                cleanedIsotopePairs.Add(pair.IsotopePairs[0]);
                            }
                        }
                        List<IsotopePair> sortedCleanedIsotopePairs = cleanedIsotopePairs.OrderBy(isotope => isotope.Isotope).ToList();
                        pair.IsotopePairs = sortedCleanedIsotopePairs;
                        isotopicDistributionFilteredPairs.Add(pair);
                    }
                    else
                    {
                        isotopicDistributionFilteredPairs.Add(pair);
                    }
                }
            }
            allPairs[rawFile] = isotopicDistributionFilteredPairs;
        }

        public List<IsotopePair> findMedianPairs(List<IsotopePair> originalPairs)
        {
            List<IsotopePair> medianPairs = new List<IsotopePair>();
            List<IsotopePair>[] sortedPairsByChannel = new List<IsotopePair>[numChannels];

            foreach (IsotopePair currentPair in originalPairs)
            {
                if (currentPair.TotalPeakCount > 0)
                {
                    for (int c = 0; c < numChannels; c++)
                    {
                        if (currentPair.MeanNormalizedIntensities[c] > 0)
                        {
                            sortedPairsByChannel[c] = originalPairs.OrderBy(pair => pair.MeanNormalizedIntensities[c]).ToList();
                        }
                    }
                }
            }

            for (int c = 0; c < numChannels; c++)
            {
                int isotopePairCount = sortedPairsByChannel[c].Count;
                IsotopePair medianIsotopePair = sortedPairsByChannel[c].ElementAt(isotopePairCount / 2);
                if (!medianPairs.Contains(medianIsotopePair)) medianPairs.Add(medianIsotopePair);
                if (isotopePairCount % 2 == 0)
                {
                    medianIsotopePair = sortedPairsByChannel[c].ElementAt((isotopePairCount / 2) - 1);
                    if (!medianPairs.Contains(medianIsotopePair)) medianPairs.Add(medianIsotopePair);
                }
            }
            return medianPairs;
        }

        /* Assembles all the peptide's quantitative information
         */
        public void quantify()
        {
            noQuantReason = NonQuantifiableType.WrongElutionProfiles;
            List<IsotopePair> isotopePairsToQuantify = new List<IsotopePair>();
            List<IsotopePair>[] isotopesToQuantify = new List<IsotopePair>[numIsotopes];

            if (Form1.ISOTOPEQUANT)
            {
                for (int i = 0; i < numIsotopes; i++)
                {
                    isotopesToQuantify[i] = new List<IsotopePair>();
                }
            }

            // NeuCode
            if (numIsotopologues > 1)
            {
                int[,] final = new int[numClusters, 2];
                if (completePairs != null && completePairs.Count > 0)
                {
                    foreach (List<Pair> pairs in completePairs.Values)
                    {
                        foreach (Pair pair in pairs)
                        {
                            foreach (IsotopePair isotopePair in pair.IsotopePairs)
                            {
                                for (int c = 0; c < numClusters; c++)
                                {
                                    if (Form1.QUANTFILTER)
                                    {
                                        // Eliminate low-level pairs by quantitative filtering
                                        if (quantFilter(isotopePair, c, true))
                                        {
                                            //int channelIndex = c * numIsotopologues;
                                            //for (int i = channelIndex; i < channelIndex + numIsotopologues; i++)
                                            //{
                                            //    double denormalizedIntensity = isotopePair.ChannelPeaks[i].GetDenormalizedIntensity(pair.InjectionTime);
                                            //    completeTotalIntensity[i, isotopePair.Isotope] += denormalizedIntensity;
                                            //    completeTotalIntensity[i, numIsotopes] += denormalizedIntensity;
                                            //}
                                            isotopePairsToQuantify.Add(isotopePair);
                                            final[c, 1]++;
                                        }
                                    }
                                    else if (Form1.ISOTOPEQUANT)
                                    {
                                        quantFilter(isotopePair, c, true);
                                        isotopesToQuantify[isotopePair.Isotope].Add(isotopePair);
                                        isotopePairsToQuantify.Add(isotopePair);
                                        final[c, 1]++;
                                    }
                                    //else
                                    //{
                                    //    int channelIndex = c * numIsotopologues;
                                    //    bool noNullPeaks = true;
                                    //    for (int i = channelIndex; i < channelIndex + numIsotopologues; i++)
                                    //    {
                                    //        if (isotopePair.ChannelPeaks[i] == null)
                                    //        {
                                    //            noNullPeaks = false;
                                    //        }
                                    //    }

                                    //    if (noNullPeaks)
                                    //    {
                                    //        for (int i = channelIndex; i < channelIndex + numIsotopologues; i++)
                                    //        {
                                    //            double denormalizedIntensity = isotopePair.ChannelPeaks[i].GetDenormalizedIntensity(isotopePair.Parent.InjectionTime);
                                    //            completeTotalIntensity[i, isotopePair.Isotope] += denormalizedIntensity;
                                    //            completeTotalIntensity[i, numIsotopes] += denormalizedIntensity;
                                    //        }
                                    //        isotopePairsToQuantify.Add(isotopePair);
                                    //        final[c, 1]++;
                                    //    }
                                    //}
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
                            foreach (IsotopePair isotopePair in pair.IsotopePairs)
                            {
                                for (int c = 0; c < numClusters; c++)
                                {
                                    int channelIndex = c * numIsotopologues;
                                    if (Form1.QUANTFILTER)
                                    {
                                        //Use a peak's intensity if it is not noise-band capped and its intensity is greater than 1/2e of the maximum intensity
                                        if (quantFilter(isotopePair, c, false))
                                        {
                                            //for (int i = channelIndex; i < channelIndex + numIsotopologues; i++)
                                            //{
                                            //    double denormalizedIntensity = isotopePair.ChannelPeaks[i].GetDenormalizedIntensity(isotopePair.Parent.InjectionTime);
                                            //    totalIntensity[i, isotopePair.Isotope] += denormalizedIntensity;
                                            //    totalIntensity[i, numIsotopes] += denormalizedIntensity;
                                            //}
                                            isotopePairsToQuantify.Add(isotopePair);
                                            final[c, 0]++;
                                        }
                                    }
                                    else if (Form1.ISOTOPEQUANT)
                                    {
                                        quantFilter(isotopePair, c, false);
                                        isotopesToQuantify[isotopePair.Isotope].Add(isotopePair);
                                        isotopePairsToQuantify.Add(isotopePair);
                                        final[c, 0]++;
                                    }
                                    //else
                                    //{
                                    //    bool noNullPeaks = true;
                                    //    int realPeaks = 0;
                                    //    for (int i = channelIndex; i < channelIndex + numIsotopologues; i++)
                                    //    {
                                    //        if (isotopePair.ChannelPeaks[i] == null)
                                    //        {
                                    //            noNullPeaks = false;
                                    //        }
                                    //        else
                                    //        {
                                    //            realPeaks++;
                                    //        }
                                    //    }

                                    //    if (noNullPeaks)
                                    //    {
                                    //        for (int i = channelIndex; i < channelIndex + numIsotopologues; i++)
                                    //        {
                                    //            if (isotopePair.ChannelPeaks[i] != null)
                                    //            {
                                    //                double denormalizedIntensity = isotopePair.ChannelPeaks[i].GetDenormalizedIntensity(isotopePair.Parent.InjectionTime);
                                    //                totalIntensity[i, isotopePair.Isotope] += denormalizedIntensity;
                                    //                totalIntensity[i, numIsotopes] += denormalizedIntensity;
                                    //            }
                                    //        }
                                    //        isotopePairsToQuantify.Add(isotopePair);
                                    //        final[c, 0]++;
                                    //    }
                                    //}
                                }
                            }
                        }
                    }
                }

                quantifiedNoiseIncluded = new bool[numClusters];
                finalQuantified = new int[numClusters];
                for (int c = 0; c < numClusters; c++)
                {
                    quantifiedNoiseIncluded[c] = true;
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

                if (Form1.LYSINEPURITYCORRECTION)
                {
                    for (int i = 0; i < numChannels; i++)
                    {
                        for (int n = 0; n <= numIsotopes; n++)
                        {
                            completeTotalIntensity[i, n] = completeTotalIntensity[i, n] / Form1.CHANNELIMPURITIES[i];
                            totalIntensity[i, n] = totalIntensity[i, n] / Form1.CHANNELIMPURITIES[i];
                        }
                    }
                }

                // Find the median ratio pairs and use for quantification
                if (final[0, 1] >= minimumTotalPairs)
                {
                    //List<IsotopePair> assembledCompletePairs = new List<IsotopePair>();
                    //foreach (List<Pair> complete in completePairs.Values)
                    //{
                    //    foreach (Pair completePair in complete)
                    //    {
                    //        foreach (IsotopePair completeIsotopePair in completePair.IsotopePairs)
                    //        {
                    //            assembledCompletePairs.Add(completeIsotopePair);
                    //        }
                    //    }
                    //}
                    List<IsotopePair> pairsToQuantify = findMedianPairs(isotopePairsToQuantify);
                    medianCompleteTotalIntensity = new double[numChannels];

                    foreach (IsotopePair quantifiedPair in pairsToQuantify)
                    {
                        for (int c = 0; c < numChannels; c++)
                        {
                            medianCompleteTotalIntensity[c] = quantifiedPair.ChannelPeaks[c].GetDenormalizedIntensity(quantifiedPair.Parent.InjectionTime) / Form1.CHANNELIMPURITIES[c];
                        }
                    }

                    for (int c = 0; c < numClusters; c++)
                    {
                        int channelIndex = c * numIsotopologues;
                        finalQuantified[c] = final[c, 1];
                        quantifiedNoiseIncluded[c] = false;
                        noQuantReason = NonQuantifiableType.Quantified;

                        //mergeCompleteAllLists();
                        foreach (MSDataFile rawFile in completePairs.Keys)
                        {
                            allPairs[rawFile] = completePairs[rawFile];
                        }

                        // Use only complete pairs for quantification
                        for (int i = channelIndex; i < channelIndex + numIsotopologues; i++)
                        {
                            for (int s = 0; s <= numIsotopes; s++)
                            {
                                totalIntensity[i, s] = completeTotalIntensity[i, s];
                            }
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
                                heavyToLightRatioSum[n - (c + 1)] = heavyInt / lightInt;
                            }

                            lightInt = medianCompleteTotalIntensity[index1];
                            heavyInt = medianCompleteTotalIntensity[n];

                            if (lightInt > 0 && heavyInt > 0)
                            {
                                heavyToLightRatioMedian[n - (c + 1), 0] = heavyInt / lightInt;
                            }
                            index2++;
                        }
                    }
                }
                else if (final[0, 1] + final[0, 0] >= minimumTotalPairs)
                {
                    //List<IsotopePair> assembledPairs = new List<IsotopePair>();
                    //foreach (List<Pair> all in allPairs.Values)
                    //{
                    //    foreach (Pair pair in all)
                    //    {
                    //        foreach (IsotopePair isotopePair in pair.IsotopePairs)
                    //        {
                    //            assembledPairs.Add(isotopePair);
                    //        }
                    //    }
                    //}
                    //List<IsotopePair> pairsToQuantify = findMedianPairs(assembledPairs);

                    List<IsotopePair> pairsToQuantify = findMedianPairs(isotopePairsToQuantify);
                    medianTotalIntensity = new double[numChannels];

                    foreach (IsotopePair quantifiedPair in pairsToQuantify)
                    {
                        for (int c = 0; c < numChannels; c++)
                        {
                            medianTotalIntensity[c] = quantifiedPair.ChannelPeaks[c].GetDenormalizedIntensity(quantifiedPair.Parent.InjectionTime) / Form1.CHANNELIMPURITIES[c];
                        }
                    }

                    mergeCompleteAllLists();
                    for (int c = 0; c < numClusters; c++)
                    {
                        int channelIndex = c * numIsotopologues;

                        for (int i = channelIndex; i < channelIndex + numIsotopologues; i++)
                        {
                            for (int s = 0; s <= numIsotopes; s++)
                            {
                                totalIntensity[i, s] += completeTotalIntensity[i, s];
                            }
                        }

                        finalQuantified[c] = final[c, 1] + final[c, 0];
                        quantifiedNoiseIncluded[c] = true;
                        noQuantReason = NonQuantifiableType.Quantified;

                        int index1 = channelIndex;
                        int index2 = index1 + 1;
                        // Use only complete pairs for quantification
                        for (int n = index2; n < channelIndex + numIsotopologues; n++)
                        {
                            lightInt = totalIntensity[index1, numIsotopes];
                            heavyInt = totalIntensity[n, numIsotopes];

                            if (lightInt > 0 && heavyInt > 0)
                            {
                                heavyToLightRatioSum[n - (c + 1)] = heavyInt / lightInt;
                            }

                            lightInt = medianTotalIntensity[index1];
                            heavyInt = medianTotalIntensity[n];

                            if (lightInt > 0 && heavyInt > 0)
                            {
                                heavyToLightRatioMedian[n - (c + 1), 0] = heavyInt / lightInt;
                            }

                            index2++;
                        }
                    }
                }
                else
                {
                    if (countAllIsotopes[0] + countCompleteIsotopes[0] < minimumTotalPairs)
                    {
                        noQuantReason = NonQuantifiableType.NotEnoughMeasurements;
                    }
                    // Not able to quantify
                    for (int i = 0; i < numChannels; i++)
                    {
                        for (int s = 0; s <= numIsotopes; s++)
                        {
                            completeTotalIntensity[i, s] = 0;
                            totalIntensity[i, s] = 0;
                        }
                    }
                }

                // First try to quantify based on complete sets of isotopologues
                /*if (completePairs != null && completePairs.Count > 0)
                {
                    foreach (List<Pair> pairs in completePairs.Values)
                    {
                        foreach (Pair pair in pairs)
                        {
                            foreach (IsotopePair isotopePair in pair.IsotopePairs)
                            {
                                for (int c = 0; c < numClusters; c++)
                                {
                                    if (Form1.QUANTFILTER)
                                    {
                                        // Eliminate low-level pairs by quantitative filtering
                                        if (quantFilter(isotopePair, c, true))
                                        {
                                            int channelIndex = c * numIsotopologues;
                                            for (int i = channelIndex; i < channelIndex + numIsotopologues; i++)
                                            {
                                                double denormalizedIntensity = isotopePair.ChannelPeaks[i].GetDenormalizedIntensity(pair.InjectionTime);
                                                completeTotalIntensity[i, isotopePair.Isotope] += denormalizedIntensity;
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
                                            if (isotopePair.ChannelPeaks[i] == null)
                                            {
                                                noNullPeaks = false;
                                            }
                                        }

                                        if (noNullPeaks)
                                        {
                                            for (int i = channelIndex; i < channelIndex + numIsotopologues; i++)
                                            {
                                                double denormalizedIntensity = isotopePair.ChannelPeaks[i].GetDenormalizedIntensity(isotopePair.Parent.InjectionTime);
                                                completeTotalIntensity[i, isotopePair.Isotope] += denormalizedIntensity;
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
                            foreach (IsotopePair isotopePair in pair.IsotopePairs)
                            {
                                for (int c = 0; c < numClusters; c++)
                                {
                                    int channelIndex = c * numIsotopologues;
                                    if (Form1.QUANTFILTER)
                                    {
                                        //Use a peak's intensity if it is not noise-band capped and its intensity is greater than 1/2e of the maximum intensity

                                        if (quantFilter(isotopePair, c, false))
                                        {
                                            for (int i = channelIndex; i < channelIndex + numIsotopologues; i++)
                                            {
                                                double denormalizedIntensity = isotopePair.ChannelPeaks[i].GetDenormalizedIntensity(isotopePair.Parent.InjectionTime);
                                                totalIntensity[i, isotopePair.Isotope] += denormalizedIntensity;
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
                                            if (isotopePair.ChannelPeaks[i] == null)
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
                                                if (isotopePair.ChannelPeaks[i] != null)
                                                {
                                                    double denormalizedIntensity = isotopePair.ChannelPeaks[i].GetDenormalizedIntensity(isotopePair.Parent.InjectionTime);
                                                    totalIntensity[i, isotopePair.Isotope] += denormalizedIntensity;
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
                }*/

                // Purity correct experiments
                /*if (Form1.LYSINEPURITYCORRECTION)
                {
                    for (int i = 0; i < numChannels; i++)
                    {
                        for (int n = 0; n <= numIsotopes; n++)
                        {
                            totalIntensity[i, n] = totalIntensity[i, n] / Form1.CHANNELIMPURITIES[i];
                            completeTotalIntensity[i, n] = completeTotalIntensity[i, n] / Form1.CHANNELIMPURITIES[i];
                        }
                    }
                } */

                /*if (coalescenceDetected)
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

                            if (final[0, 1] > 0) mergeCompleteAllLists();
                            else
                            {
                                foreach (MSDataFile rawFile in completePairs.Keys)
                                {
                                    allPairs[rawFile] = completePairs[rawFile];
                                }
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
                                lightInt = totalIntensity[index1, numIsotopes];
                                heavyInt = totalIntensity[n, numIsotopes];

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

                            if (final[0, 1] > 0) mergeCompleteAllLists();

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
                }*/
            }

            // Traditional SILAC

            else
            {
                int[,] final = new int[1, 2];
                if (completePairs != null && completePairs.Count > 0)
                {
                    foreach (List<Pair> pairs in completePairs.Values)
                    {
                        foreach (Pair pair in pairs)
                        {
                            foreach (IsotopePair isotopePair in pair.IsotopePairs)
                            {
                                    if (Form1.QUANTFILTER)
                                    {
                                        // Eliminate low-level pairs by quantitative filtering
                                        if (quantFilter(isotopePair, 0, true))
                                        {
                                            //for (int i = 0; i < numChannels; i++)
                                            //{
                                            //    double denormalizedIntensity = isotopePair.ChannelPeaks[i].GetDenormalizedIntensity(pair.InjectionTime);
                                            //    completeTotalIntensity[i, isotopePair.Isotope] += denormalizedIntensity;
                                            //    completeTotalIntensity[i, numIsotopes] += denormalizedIntensity;
                                            //}
                                            isotopePairsToQuantify.Add(isotopePair);
                                            final[0, 1]++;
                                        }
                                    }
                                    else if (Form1.ISOTOPEQUANT)
                                    {
                                        quantFilter(isotopePair, 0, true);
                                        isotopesToQuantify[isotopePair.Isotope].Add(isotopePair);
                                        isotopePairsToQuantify.Add(isotopePair);
                                        final[0, 1]++;
                                    }
                                    //else
                                    //{
                                    //    bool noNullPeaks = true;
                                    //    for (int i = 0; i < numChannels; i++)
                                    //    {
                                    //        if (isotopePair.ChannelPeaks[i] == null)
                                    //        {
                                    //            noNullPeaks = false;
                                    //        }
                                    //    }

                                    //    if (noNullPeaks)
                                    //    {
                                    //        for (int i = 0; i < numChannels; i++)
                                    //        {
                                    //            double denormalizedIntensity = isotopePair.ChannelPeaks[i].GetDenormalizedIntensity(isotopePair.Parent.InjectionTime);
                                    //            completeTotalIntensity[i, isotopePair.Isotope] += denormalizedIntensity;
                                    //            completeTotalIntensity[i, numIsotopes] += denormalizedIntensity;
                                    //        }
                                    //        isotopePairsToQuantify.Add(isotopePair);
                                    //        final[0, 1]++;
                                    //    }
                                    //}
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
                            foreach (IsotopePair isotopePair in pair.IsotopePairs)
                            {
                                if (Form1.QUANTFILTER)
                                {
                                    //Use a peak's intensity if it is not noise-band capped and its intensity is greater than 1/2e of the maximum intensity
                                    if (quantFilter(isotopePair, 0, false))
                                    {
                                        //for (int i = 0; i < numChannels; i++)
                                        //{
                                        //    double denormalizedIntensity = isotopePair.ChannelPeaks[i].GetDenormalizedIntensity(isotopePair.Parent.InjectionTime);
                                        //    totalIntensity[i, isotopePair.Isotope] += denormalizedIntensity;
                                        //    totalIntensity[i, numIsotopes] += denormalizedIntensity;

                                        //}
                                        isotopePairsToQuantify.Add(isotopePair);
                                        final[0, 0]++;
                                    }
                                }
                                else if (Form1.ISOTOPEQUANT)
                                {
                                    quantFilter(isotopePair, 0, false);
                                    isotopesToQuantify[isotopePair.Isotope].Add(isotopePair);
                                    isotopePairsToQuantify.Add(isotopePair);
                                    final[0, 0]++;
                                }
                                //else
                                //{
                                //    bool noNullPeaks = true;
                                //    int realPeaks = 0;
                                //    for (int i = 0; i < numChannels; i++)
                                //    {
                                //        if (isotopePair.ChannelPeaks[i] == null)
                                //        {
                                //            noNullPeaks = false;
                                //        }
                                //        else
                                //        {
                                //            realPeaks++;
                                //        }
                                //    }

                                //    if (noNullPeaks)
                                //    {
                                //        for (int i = 0; i < numChannels; i++)
                                //        {
                                //            if (isotopePair.ChannelPeaks[i] != null)
                                //            {
                                //                double denormalizedIntensity = isotopePair.ChannelPeaks[i].GetDenormalizedIntensity(isotopePair.Parent.InjectionTime);
                                //                totalIntensity[i, isotopePair.Isotope] += denormalizedIntensity;
                                //                totalIntensity[i, numIsotopes] += denormalizedIntensity;
                                //            }
                                //        }
                                //        isotopePairsToQuantify.Add(isotopePair);
                                //        final[0, 0]++;
                                //    }
                                //}
                            }
                        }
                    }
                }

                quantifiedNoiseIncluded = new bool[1];
                quantifiedNoiseIncluded[0] = true;
                finalQuantified = new int[1];

                double lightInt;
                double heavyInt;
                bool quantified = false;
                int minimumTotalPairs = 3;
                int minimumNoNBCPairs = 3;
                int minimumPostQFPairs = 3;

                if (Form1.LYSINEPURITYCORRECTION)
                {
                    for (int i = 0; i < numChannels; i++)
                    {
                        for (int n = 0; n <= numIsotopes; n++)
                        {
                            completeTotalIntensity[i, n] = completeTotalIntensity[i, n] / Form1.CHANNELIMPURITIES[i];
                            totalIntensity[i, n] = totalIntensity[i, n] / Form1.CHANNELIMPURITIES[i];
                        }
                    }
                }

                // Find the median ratio pairs and use for quantification
                if (final[0, 1] >= minimumTotalPairs)
                {
                    //mergeCompleteAllLists();
                    //List<IsotopePair> assembledCompletePairs = new List<IsotopePair>();
                    //foreach (List<Pair> complete in completePairs.Values)
                    //{
                    //    foreach (Pair completePair in complete)
                    //    {
                    //        foreach (IsotopePair completeIsotopePair in completePair.IsotopePairs)
                    //        {
                    //            assembledCompletePairs.Add(completeIsotopePair);
                    //        }
                    //    }
                    //}
                    List<IsotopePair> pairsToQuantify = findMedianPairs(isotopePairsToQuantify);
                    medianCompleteTotalIntensity = new double[numChannels];

                    foreach (IsotopePair quantifiedPair in pairsToQuantify)
                    {
                        for (int c = 0; c < numChannels; c++)
                        {
                            medianCompleteTotalIntensity[c] += quantifiedPair.ChannelPeaks[c].GetDenormalizedIntensity(quantifiedPair.Parent.InjectionTime);
                        }
                    }

                    finalQuantified[0] = final[0, 1];
                    quantifiedNoiseIncluded[0] = false;
                    noQuantReason = NonQuantifiableType.Quantified;

                    foreach (MSDataFile rawFile in completePairs.Keys)
                    {
                        allPairs[rawFile] = completePairs[rawFile];
                    }

                    // Use only complete pairs for quantification
                    for (int i = 0; i < numChannels; i++)
                    {
                        for (int s = 0; s <= numIsotopes; s++)
                        {
                            totalIntensity[i, s] = completeTotalIntensity[i, s];
                        }
                    }

                    int index1 = 0;
                    int index2 = index1 + 1;
                    // Use only complete pairs for quantification
                    for (int n = index2; n < numChannels; n++)
                    {
                        lightInt = totalIntensity[index1, numIsotopes];
                        heavyInt = totalIntensity[n, numIsotopes];

                        if (lightInt > 0 && heavyInt > 0)
                        {
                            heavyToLightRatioSum[n - 1] = heavyInt / lightInt;
                        }
                        index2++;
                    }
                }
                else if (final[0, 1] + final[0, 0] >= minimumTotalPairs)
                {
                    mergeCompleteAllLists();

                    //List<IsotopePair> assembledPairs = new List<IsotopePair>();
                    //foreach (List<Pair> all in allPairs.Values)
                    //{
                    //    foreach (Pair pair in all)
                    //    {
                    //        foreach (IsotopePair isotopePair in pair.IsotopePairs)
                    //        {
                    //            assembledPairs.Add(isotopePair);
                    //        }
                    //    }
                    //}
                    //List<IsotopePair> pairsToQuantify = findMedianPairs(assembledPairs);

                    List<IsotopePair> pairsToQuantify = findMedianPairs(isotopePairsToQuantify);
                    medianTotalIntensity = new double[numChannels];

                    foreach (IsotopePair quantifiedPair in pairsToQuantify)
                    {
                        for (int c = 0; c < numChannels; c++)
                        {
                            medianTotalIntensity[c] += quantifiedPair.ChannelPeaks[c].GetDenormalizedIntensity(quantifiedPair.Parent.InjectionTime);
                        }
                    }

                    for (int i = 0; i < numChannels; i++)
                    {
                        for (int s = 0; s <= numIsotopes; s++)
                        {
                            totalIntensity[i, s] += completeTotalIntensity[i, s];
                        }
                    }

                    finalQuantified[0] = final[0, 1] + final[0, 0];
                    quantifiedNoiseIncluded[0] = true;
                    noQuantReason = NonQuantifiableType.Quantified;

                    int index1 = 0;
                    int index2 = index1 + 1;
                    // Use only complete pairs for quantification
                    for (int n = index2; n < numChannels; n++)
                    {
                        lightInt = totalIntensity[index1, numIsotopes];
                        heavyInt = totalIntensity[n, numIsotopes];

                        if (lightInt > 0 && heavyInt > 0)
                        {
                            heavyToLightRatioSum[n - 1] = heavyInt / lightInt;
                        }
                        index2++;
                    }
                }
                else
                {
                    if (countAllIsotopes[0] + countCompleteIsotopes[0] < minimumTotalPairs)
                    {
                        noQuantReason = NonQuantifiableType.NotEnoughMeasurements;
                    }
                    // Not able to quantify
                    for (int i = 0; i < numChannels; i++)
                    {
                        for (int s = 0; s <= numIsotopes; s++)
                        {
                            completeTotalIntensity[i, s] = 0;
                            totalIntensity[i, s] = 0;
                        }
                    }
                }
            }
            //else
            //{
            //    int[,] final = new int[1, 2];
            //    quantifiedNoiseIncluded = new bool[1];
            //    finalQuantified = new int[1];

            //    double lightInt;
            //    double heavyInt;
            //    bool quantified = false;
            //    int minimumTotalPairs;
            //    int minimumNoNBCPairs;
            //    int minimumPostQFPairs;

            //    if (Form1.MULTIINJECT)
            //    {
            //        minimumTotalPairs = 1;
            //        minimumNoNBCPairs = 1;
            //        minimumPostQFPairs = 1;
            //    }
            //    else
            //    {
            //        minimumTotalPairs = 3;
            //        minimumNoNBCPairs = 3;
            //        minimumPostQFPairs = 3;
            //    }

            //    // First try to quantify based on complete sets of isotopologues
            //    if (completePairs != null && completePairs.Count > 0)
            //    {
            //        foreach (List<Pair> pairs in completePairs.Values)
            //        {
            //            foreach (Pair pair in pairs)
            //            {
            //                foreach (IsotopePair isotopePair in pair.IsotopePairs)
            //                {
            //                    if (Form1.QUANTFILTER)
            //                    {
            //                        // Eliminate low-level pairs by quantitative filtering
            //                        if (quantFilter(isotopePair, 0, true))
            //                        {
            //                            for (int i = 0; i < numChannels; i++)
            //                            {
            //                                double denormalizedIntensity = isotopePair.ChannelPeaks[i].GetDenormalizedIntensity(isotopePair.Parent.InjectionTime);
            //                                completeTotalIntensity[i, isotopePair.Isotope] += denormalizedIntensity;
            //                                completeTotalIntensity[i, numIsotopes] += denormalizedIntensity;
            //                            }
            //                            final[0, 1]++;
            //                        }
            //                    }
            //                    else
            //                    {
            //                        bool noNullPeaks = true;
            //                        for (int i = 0; i < numChannels; i++)
            //                        {
            //                            if (isotopePair.ChannelPeaks[i] == null)
            //                            {
            //                                noNullPeaks = false;
            //                            }
            //                        }

            //                        if (noNullPeaks)
            //                        {
            //                            for (int i = 0; i < numChannels; i++)
            //                            {
            //                                double denormalizedIntensity = isotopePair.ChannelPeaks[i].GetDenormalizedIntensity(isotopePair.Parent.InjectionTime);
            //                                completeTotalIntensity[i, isotopePair.Isotope] += denormalizedIntensity;
            //                                completeTotalIntensity[i, numIsotopes] += denormalizedIntensity;
            //                            }
            //                            final[0, 1]++;
            //                        }
            //                    }
            //                } // End isotope loop
            //            } // End pair loop
            //        } // End pair list loop
            //    }
            //    if (allPairs != null && allPairs.Count > 0)
            //    {
            //        foreach (List<Pair> pairs in allPairs.Values)
            //        {
            //            foreach (Pair pair in pairs)
            //            {
            //                foreach (IsotopePair isotopePair in pair.IsotopePairs)
            //                {
            //                    if (Form1.QUANTFILTER)
            //                    {
            //                        //Use a peak's intensity if it is not noise-band capped and its intensity is greater than 1/2e of the maximum intensity
            //                        if (quantFilter(isotopePair, 0, false))
            //                        {
            //                            for (int i = 0; i < numChannels; i++)
            //                            {
            //                                double denormalizedIntensity = isotopePair.ChannelPeaks[i].GetDenormalizedIntensity(isotopePair.Parent.InjectionTime);
            //                                totalIntensity[i, isotopePair.Isotope] += denormalizedIntensity;
            //                                totalIntensity[i, numIsotopes] += denormalizedIntensity;

            //                            }
            //                            final[0, 0]++;
            //                        }
            //                    }
            //                    else
            //                    {
            //                        bool noNullPeaks = true;
            //                        int realPeaks = 0;
            //                        for (int i = 0; i < numChannels; i++)
            //                        {
            //                            if (isotopePair.ChannelPeaks[i] == null)
            //                            {
            //                                noNullPeaks = false;
            //                            }
            //                            else
            //                            {
            //                                realPeaks++;
            //                            }
            //                        }

            //                        if (noNullPeaks)
            //                        {
            //                            for (int i = 0; i < numChannels; i++)
            //                            {
            //                                if (isotopePair.ChannelPeaks[i] != null)
            //                                {
            //                                    double denormalizedIntensity = isotopePair.ChannelPeaks[i].GetDenormalizedIntensity(isotopePair.Parent.InjectionTime);
            //                                    totalIntensity[i, isotopePair.Isotope] += denormalizedIntensity;
            //                                    totalIntensity[i, numIsotopes] += denormalizedIntensity;
            //                                }
            //                            }
            //                            final[0, 0]++;
            //                        }
            //                    }
            //                }
            //            }
            //        }
            //    }

            //    // Purity correct experiments
            //    if (Form1.LYSINEPURITYCORRECTION)
            //    {
            //        for (int i = 0; i < numChannels; i++)
            //        {
            //            for (int n = 0; n <= numIsotopes; n++)
            //            {
            //                totalIntensity[i, n] = totalIntensity[i, n] / Form1.CHANNELIMPURITIES[i];
            //                completeTotalIntensity[i, n] = completeTotalIntensity[i, n] / Form1.CHANNELIMPURITIES[i];
            //            }
            //        }
            //    }

            //    if (final[0, 1] >= minimumPostQFPairs)
            //    {
            //        finalQuantified[0] = final[0, 1];
            //        quantifiedNoiseIncluded[0] = false;
            //        noQuantReason = NonQuantifiableType.Quantified;

            //        if (final[0, 1] > 0) mergeCompleteAllLists();
            //        else
            //        {
            //            foreach (MSDataFile rawFile in completePairs.Keys)
            //            {
            //                allPairs[rawFile] = completePairs[rawFile];
            //            }
            //        }

            //        // Use only complete pairs for quantification
            //        for (int i = 0; i < numChannels; i++)
            //        {
            //            totalIntensity[i, numIsotopes] += completeTotalIntensity[i, numIsotopes];
            //        }
                    
            //        int index1 = 0;
            //        int index2 = index1 + 1;
                    
            //        for (int i = index2; i < numChannels; i++)
            //        {
            //            lightInt = totalIntensity[index1, numIsotopes];
            //            heavyInt = totalIntensity[i, numIsotopes];

            //            if (lightInt > 0 && heavyInt > 0)
            //            {
            //                heavyToLightRatioSum[i - 1, 0] = heavyInt / lightInt;
            //            }
            //        }
            //    }
            //    else if (final[0, 0] + final[0, 1] >= minimumPostQFPairs)
            //    {
            //        finalQuantified[0] = final[0, 0] + final[0, 1];
            //        quantifiedNoiseIncluded[0] = true;
            //        noQuantReason = NonQuantifiableType.Quantified;

            //        if (final[0, 1] > 0) mergeCompleteAllLists();
                    
            //        int index1 = 0;
            //        int index2 = index1 + 1;

            //        // Use all pairs for quantification
            //        for (int i = 0; i < numChannels; i++)
            //        {
            //            totalIntensity[i, numIsotopes] += completeTotalIntensity[i, numIsotopes];
            //        }
                    
            //        for (int i = index2; i < numChannels; i++)
            //        {
            //            lightInt = totalIntensity[index1, numIsotopes];
            //            heavyInt = totalIntensity[i, numIsotopes];

            //            if (lightInt > 0 && heavyInt > 0)
            //            {
            //                heavyToLightRatioSum[i - 1, 0] = heavyInt / lightInt;
            //            }
            //        }
            //    }
            //    else
            //    {
            //        if (countAllIsotopes[0] + countCompleteIsotopes[0] >= minimumPostQFPairs)
            //        {
            //            noQuantReason = NonQuantifiableType.WrongElutionProfiles;
            //        }
            //        else if (countAllIsotopes[0] + countCompleteIsotopes[0] > 0)
            //        {
            //            noQuantReason = NonQuantifiableType.NotEnoughMeasurements;
            //        }

            //        // Not able to quantify
            //        for (int i = 0; i < numChannels; i++)
            //        {
            //            completeTotalIntensity[i, numIsotopes] = 0;
            //            totalIntensity[i, numIsotopes] = 0;
            //        }
            //    }

            //    //int[,] final = new int[1, 2];
            //    //quantifiedNoiseIncluded = new bool[1];
            //    //finalQuantified = new int[1];

            //    //if (completePairs != null && completePairs.Count > 0)
            //    //{
            //    //    foreach (List<Pair> pairs in completePairs.Values)
            //    //    {
            //    //        foreach (Pair pair in pairs)
            //    //        {
            //    //            for (int j = 0; j < numIsotopes; j++)
            //    //            {
            //    //                if (Form1.QUANTFILTER)
            //    //                {
            //    //                    // Eliminate low-level pairs by quantitative filtering
            //    //                    if (quantFilter(pair, j, 0, true))
            //    //                    {
            //    //                        for (int i = 0; i < numChannels; i++)
            //    //                        {
            //    //                            double denormalizedIntensity = pair.peaks[i, j].GetDenormalizedIntensity(pair.injectionTime);
            //    //                            completeTotalIntensity[i, j] += denormalizedIntensity;
            //    //                            completeTotalIntensity[i, numIsotopes] += denormalizedIntensity;
            //    //                        }
            //    //                        final[0, 1]++;
            //    //                    }
            //    //                }
            //    //                else
            //    //                {
            //    //                    bool noNullPeaks = true;
            //    //                    for (int i = 0; i < numChannels; i++)
            //    //                    {
            //    //                        if (pair.peaks[i, j] == null)
            //    //                        {
            //    //                            noNullPeaks = false;
            //    //                        }
            //    //                    }

            //    //                    if (noNullPeaks)
            //    //                    {
            //    //                        for (int i = 0; i < numChannels; i++)
            //    //                        {
            //    //                            double denormalizedIntensity = pair.peaks[i, j].GetDenormalizedIntensity(pair.injectionTime);
            //    //                            completeTotalIntensity[i, j] += denormalizedIntensity;
            //    //                            completeTotalIntensity[i, numIsotopes] += denormalizedIntensity;
            //    //                        }
            //    //                        final[0, 1]++;
            //    //                    }
            //    //                }
            //    //            } // End isotope loop
            //    //        } // End pair loop
            //    //    } // End pair list loop
            //    //}
            //    //if (allPairs != null && allPairs.Count > 0)
            //    //{
            //    //    foreach (List<Pair> pairs in allPairs.Values)
            //    //    {
            //    //        foreach (Pair pair in pairs)
            //    //        {
            //    //            for (int j = 0; j < numIsotopes; j++)
            //    //            {
            //    //                if (Form1.QUANTFILTER)
            //    //                {
            //    //                    //Use a peak's intensity if it is not noise-band capped and its intensity is greater than 1/2e of the maximum intensity

            //    //                    if (quantFilter(pair, j, 0, false))
            //    //                    {
            //    //                        for (int i = 0; i < numChannels; i++)
            //    //                        {
            //    //                            double denormalizedIntensity = pair.peaks[i, j].GetDenormalizedIntensity(pair.injectionTime);
            //    //                            totalIntensity[i, j] += denormalizedIntensity;
            //    //                            totalIntensity[i, numIsotopes] += denormalizedIntensity;

            //    //                        }
            //    //                        final[0, 0]++;
            //    //                    }
            //    //                }
            //    //                else
            //    //                {
            //    //                    bool noNullPeaks = true;
            //    //                    int realPeaks = 0;
            //    //                    for (int i = 0; i < numChannels; i++)
            //    //                    {
            //    //                        if (pair.peaks[i, j] == null)
            //    //                        {
            //    //                            noNullPeaks = false;
            //    //                        }
            //    //                        else
            //    //                        {
            //    //                            realPeaks++;
            //    //                        }
            //    //                    }

            //    //                    if (noNullPeaks)
            //    //                    {
            //    //                        for (int i = 0; i < numChannels; i++)
            //    //                        {
            //    //                            if (pair.peaks[i, j] != null)
            //    //                            {
            //    //                                double denormalizedIntensity = pair.peaks[i, j].GetDenormalizedIntensity(pair.injectionTime);
            //    //                                totalIntensity[i, j] += denormalizedIntensity;
            //    //                                totalIntensity[i, numIsotopes] += denormalizedIntensity;
            //    //                            }
            //    //                        }
            //    //                        final[0, 0]++;
            //    //                    }
            //    //                }
            //    //            }
            //    //        }
            //    //    }
            //    //}

            //    //double lightInt;
            //    //double heavyInt;
            //    //bool quantified = false;
            //    //int minimumTotalPairs = 3;
            //    //int minimumNoNBCPairs = 3;
            //    //int minimumPostQFPairs = 3;

            //    //if (final[0, 1] >= minimumPostQFPairs)
            //    //{
            //    //    finalQuantified[0] = final[0, 1];
            //    //    quantifiedNoiseIncluded[0] = false;

            //    //    for (int i = 1; i < numChannels; i++)
            //    //    {
            //    //        if (Form1.NOISEBANDCAP)
            //    //        {
            //    //            totalIntensity = completeTotalIntensity;
            //    //        }
            //    //        // Use only complete pairs for quantification
            //    //        int index1 = 0;
            //    //        int index2 = i;

            //    //        lightInt = completeTotalIntensity[index1, numIsotopes];
            //    //        heavyInt = completeTotalIntensity[index2, numIsotopes];

            //    //        if (lightInt > 0 && heavyInt > 0)
            //    //        {
            //    //            heavyToLightRatioSum[i - 1, 0] = heavyInt / lightInt;
            //    //        }
            //    //    }
            //    //}
            //    //else if (final[0, 0] >= minimumPostQFPairs)
            //    //{
            //    //    finalQuantified[0] = final[0, 0];
            //    //    quantifiedNoiseIncluded[0] = true;

            //    //    for (int i = 1; i < numChannels; i++)
            //    //    {
            //    //        int index1 = 0;
            //    //        int index2 = i;

            //    //        lightInt = totalIntensity[index1, numIsotopes];
            //    //        heavyInt = totalIntensity[index2, numIsotopes];

            //    //        if (lightInt > 0 && heavyInt > 0)
            //    //        {
            //    //            heavyToLightRatioSum[i - 1, 0] = heavyInt / lightInt;
            //    //        }
            //    //    }
            //    //}
            //    //else
            //    //{
            //    //    // Not able to quantify
            //    //    for (int i = 0; i < numChannels; i++)
            //    //    {
            //    //        completeTotalIntensity[i, numIsotopes] = 0;
            //    //        totalIntensity[i, numIsotopes] = 0;
            //    //    }
            //    //}
            //}  
        }

        public void mergeCompleteAllLists()
        {
            List<Pair> mergedList;
            Pair pairToUpdate = null;
            int scanNumber;
            foreach (MSDataFile rawFile in completePairs.Keys)
            {
                mergedList = allPairs[rawFile];
                List<Pair> completePairList = completePairs[rawFile];

                foreach (Pair pairToAdd in completePairList)
                {
                    scanNumber = pairToAdd.ScanNumber;
                    bool found = false;
                    foreach (Pair currentPair in mergedList)
                    {
                        if (currentPair.ScanNumber == scanNumber)
                        {
                            found = true;
                            pairToUpdate = currentPair;
                        }
                    }

                    if (numIsotopologues > 1 && numClusters > 1)
                    {
                        for (int i = 0; i < numClusters; i++)
                        {
                            if (!found && pairToAdd.TotalPeakCount > 0) mergedList.Add(pairToAdd);
                            else
                            {
                                foreach (IsotopePair isotopePairToAdd in pairToAdd.IsotopePairs)
                                {
                                    if (isotopePairToAdd.CompleteByCluster[i])
                                    {
                                        bool isotopeFound = false;
                                        IsotopePair isotopePairToUpdate = null;
                                        foreach (IsotopePair currentIsotopePair in pairToUpdate.IsotopePairs)
                                        {
                                            if (currentIsotopePair.Isotope == isotopePairToAdd.Isotope)
                                            {
                                                isotopeFound = true;
                                                isotopePairToUpdate = currentIsotopePair;
                                            }
                                        }

                                        if (!found) pairToUpdate.IsotopePairs.Add(isotopePairToAdd);
                                        else
                                        {
                                            int channelIndex = i * numIsotopologues;
                                            for (int m = channelIndex; m < channelIndex + numIsotopologues; m++)
                                            {
                                                isotopePairToUpdate.ChannelPeaks[m] = isotopePairToAdd.ChannelPeaks[m];
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    else
                    {
                        if (!found && pairToAdd.TotalPeakCount > 0) mergedList.Add(pairToAdd);
                        else
                        {
                            foreach (IsotopePair isotopePairToAdd in pairToAdd.IsotopePairs)
                            {
                                if (isotopePairToAdd.Complete) pairToUpdate.IsotopePairs.Add(isotopePairToAdd);
                            }
                        }        
                    }

                }
                allPairs[rawFile] = mergedList;
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
