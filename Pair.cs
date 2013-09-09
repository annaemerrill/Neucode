using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using CSMSL;
using CSMSL.IO;
using CSMSL.Spectral;
using CSMSL.Proteomics;
using CSMSL.Chemistry;
using CSMSL.IO.Thermo;
using CSMSL.Proteomics;


namespace Coon.NeuQuant
{
    class Pair
    {
        public PeptideID parent; // The peptide the pair belongs to
        public int scanNumber; // The scan number the pair is associated with
        public MSDataFile rawFile; // The raw file the pair is associated with
        public int totalPeakCount
        {
            get
            {
                int count = 0;
                int channels = parent.numChannels;
                int isotopes = parent.numIsotopes;
                int clusters = parent.numClusters;
                if (peakCount != null)
                {
                    if (parent.numIsotopologues > 1)
                    {
                        for (int c = 0; c < clusters; c++)
                        {
                            for (int i = 0; i < isotopes; i++)
                            {
                                count += peakCount[c, i];
                            }
                        }
                    }
                    else
                    {
                        for (int i = 0; i < isotopes; i++)
                        {
                            count += peakCount[0, i];
                        }
                    }
                }
                return count;
            }
        } // The total number of peaks associated with this pair
        public ILabeledPeak[,] peaks; // Keeps track of peaks for quantitation (rows = # channels; columns = # isotopes)
        public bool[,] complete // Keeps track of whether or not a cluster (i.e., set of isotopologues) of peaks is complete (rows = # clusters; columns = # isotopes)
        {
            get
            {
                int channels = parent.numChannels;
                int isotopes = parent.numIsotopes;
                int clusters = parent.numClusters;
                int isotopologues = parent.numIsotopologues;
                int channelIndex;
                bool[,] complete;
                if (isotopologues > 1)
                {
                    complete = new bool[clusters, isotopes];
                    if (peaks != null)
                    {
                        for (int c = 0; c < clusters; c++)
                        {
                            for (int j = 0; j < isotopes; j++)
                            {
                                complete[c, j] = true;
                                channelIndex = c * isotopologues;
                                for (int i = channelIndex; i < channelIndex + isotopologues; i++)
                                {
                                    // A cluster is incomplete if any of its peaks are null or noise-band capped
                                    if (peaks[i, j] == null || peaks[i, j].X == 0)
                                    {
                                        complete[c, j] = false;
                                    }
                                }
                            }
                        }
                    }
                    return complete;
                }
                else
                {
                    complete = new bool[1, isotopes];
                    if (peaks != null)
                    {
                        for (int j = 0; j < isotopes; j++)
                        {
                            complete[0,j] = true;
                            for (int i = 0; i < channels; i++)
                            {
                                if (peaks[i, j] == null || peaks[i, j].X == 0)
                                {
                                    complete[0, j] = false;
                                }
                            }
                        }
                    }
                    return complete;
                }
            }
        }
        public double[,] maxClusterIntensity // Keeps track of each cluster's maximum intensity (rows = # clusters; columns = # isotopes)
        {
            get
            {
                int channels = parent.numChannels;
                int clusters = parent.numClusters;
                int isotopes = parent.numIsotopes;
                int isotopologues = parent.numIsotopologues;
                int channelIndex;
                double[,] max = new double[clusters, isotopes];
                {
                    if (peaks != null)
                    {
                        for (int j = 0; j < isotopes; j++)
                        {
                            for (int c = 0; c < clusters; c++)
                            {
                                max[c,j] = 0;
                                channelIndex = c * isotopologues;
                                for (int i = channelIndex; i < channelIndex + isotopologues; i++)
                                {
                                    if (peaks[i, j] != null && peaks[i, j].GetDenormalizedIntensity(injectionTime) > max[c,j])
                                    {
                                        max[c,j] = peaks[i, j].GetDenormalizedIntensity(injectionTime);
                                    }
                                }
                            }
                        }
                    }
                }
                return max;
            }
        }
        public int[,] peakCount // Keeps track of the number of peaks in each cluster (rows = # clusters; columns = # isotopes)
        {
            get
            {
                int channels = parent.numChannels;
                int isotopes = parent.numIsotopes;
                int clusters = parent.numClusters;
                int isotopologues = parent.numIsotopologues;
                int[,] count = null;
                if (peaks != null)
                {
                    if (isotopologues > 1)
                    {
                        count = new int[clusters, isotopes];
                        for (int j = 0; j < isotopes; j++)
                        {
                            for (int c = 0; c < clusters; c++)
                            {
                                int channelIndex = c * isotopologues;
                                for (int i = channelIndex; i < channelIndex + isotopologues; i++)
                                {
                                    if (peaks[i, j] != null)
                                    {
                                        count[c, j]++;
                                    }
                                }
                            }
                        }
                        return count;
                    }
                    else
                    {
                        count = new int[1, isotopes];
                        for (int j = 0; j < isotopes; j++)
                        {
                            for (int i = 0; i < channels; i++)
                            {
                                if (peaks[i, j] != null)
                                {
                                    count[0, j]++;
                                }
                            }
                        }
                        return count;
                    }
                }
                return null;
            }
        }
        public int[,] peakPattern
        {
            get
            {
                int channels = parent.numChannels;
                int isotopes = parent.numIsotopes;
                int clusters = parent.numClusters;
                int isotopologues = parent.numIsotopologues;
                int[,] count = null;
                if (peaks != null)
                {
                    // NeuCode quantification
                    if (isotopologues > 1)
                    {
                        double minSignal = Form1.MINIMUMSN * 1.0;
                        count = new int[clusters, isotopes];
                        for (int j = 0; j < isotopes; j++)
                        {
                            for (int c = 0; c < clusters; c++)
                            {
                                if (peakCount[c, j] == isotopologues || peakCount[c, j] < parent.peaksNeeded) count[c, j] = 0;
                                else
                                {
                                    if (isotopologues == 2)
                                    {
                                        // 1 peak pattern
                                        if (peakCount[c, j] == 1)
                                        {
                                            if (peaks[c + 1, j] == null)
                                            {
                                                // Peak 1
                                                if (peaks[c, j].GetSignalToNoise() > minSignal) count[c, j] = 0;
                                                else count[c, j] = -1;
                                            }
                                            else if (peaks[c, j] == null)
                                            {
                                                // Peak 2
                                                if (peaks[c + 1, j].GetSignalToNoise() > minSignal) count[c, j] = 1;
                                                else count[c, j] = -1;
                                            }
                                        }
                                    }
                                    else if (isotopologues == 3)
                                    {
                                        // 1 peak pattern
                                        if (peakCount[c, j] == 1)
                                        {
                                            if (peaks[c + 2, j] == null)
                                            {
                                                // Peak 1
                                                if (peaks[c + 1, j] == null && peaks[c, j].GetSignalToNoise() > minSignal) count[c, j] = 0;
                                                // Peak 2
                                                else if (peaks[c, j] == null && peaks[c + 1, j].GetSignalToNoise() > minSignal) count[c, j] = 1;
                                                else count[c, j] = -1;
                                            }
                                            // Peak 3
                                            else if (peaks[c + 2, j].GetSignalToNoise() > minSignal) count[c, j] = 2;
                                            else count[c, j] = -1;
                                        }
                                        // 2 peak pattern
                                        else if (peakCount[c, j] == 2)
                                        {
                                            // Peaks 1 & 2
                                            if (peaks[c + 2, j] == null && peaks[c, j].GetSignalToNoise() > minSignal && peaks[c + 1, j].GetSignalToNoise() > minSignal) count[c, j] = 0;
                                            // Peaks 1 & 3
                                            else if (peaks[c + 1, j] == null && peaks[c, j].GetSignalToNoise() > minSignal && peaks[c + 2, j].GetSignalToNoise() > minSignal) count[c, j] = 1;
                                            // Peaks 2 & 3
                                            else if (peaks[c, j] == null && peaks[c + 1, j].GetSignalToNoise() > minSignal && peaks[c + 2, j].GetSignalToNoise() > minSignal) count[c, j] = 2;
                                            else count[c, j] = -1;
                                        }
                                    }
                                    else if (isotopologues == 4)
                                    {
                                        // 2 peak pattern
                                        if (peakCount[c, j] == 2)
                                        {
                                            if (peaks[c + 3, j] == null)
                                            {
                                                // Peaks 1 & 2
                                                if (peaks[c + 2, j] == null && peaks[c, j].GetSignalToNoise() > minSignal && peaks[c + 1, j].GetSignalToNoise() > minSignal) count[c, j] = 0;
                                                // Peaks 1 & 3
                                                else if (peaks[c + 1, j] == null && peaks[c, j].GetSignalToNoise() > minSignal && peaks[c + 2, j].GetSignalToNoise() > minSignal) count[c, j] = 1;
                                                // Peaks 2 & 3
                                                else if (peaks[c, j] == null && peaks[c + 1, j].GetSignalToNoise() > minSignal && peaks[c + 2, j].GetSignalToNoise() > minSignal) count[c, j] = 3;
                                                else count[c, j] = -1;
                                            }
                                            else if (peaks[c + 2, j] == null)
                                            {
                                                // Peaks 1 & 4
                                                if (peaks[c + 1, j] == null && peaks[c, j].GetSignalToNoise() > minSignal && peaks[c + 3, j].GetSignalToNoise() > minSignal) count[c, j] = 2;
                                                // Peaks 2 & 4
                                                else if (peaks[c, j] == null && peaks[c + 1, j].GetSignalToNoise() > minSignal && peaks[c + 3, j].GetSignalToNoise() > minSignal) count[c, j] = 4;
                                                else count[c, j] = -1;
                                            }
                                            // Peaks 3 & 4
                                            else if (peaks[c + 1, j] == null && peaks[c + 2, j].GetSignalToNoise() > minSignal && peaks[c + 3, j].GetSignalToNoise() > minSignal) count[c, j] = 5;
                                        }
                                        // 3 peak pattern
                                        else if (peakCount[c, j] == 3)
                                        {
                                            // Peaks 1 & 2 & 3
                                            if (peaks[c + 3, j] == null && peaks[c, j].GetSignalToNoise() > minSignal && peaks[c + 1, j].GetSignalToNoise() > minSignal && peaks[c + 2, j].GetSignalToNoise() > minSignal) count[c, j] = 0;
                                            // Peaks 1 & 2 & 4
                                            else if (peaks[c + 2, j] == null && peaks[c, j].GetSignalToNoise() > minSignal && peaks[c + 1, j].GetSignalToNoise() > minSignal && peaks[c + 3, j].GetSignalToNoise() > minSignal) count[c, j] = 1;
                                            // Peaks 1 & 3 & 4
                                            else if (peaks[c + 1, j] == null && peaks[c, j].GetSignalToNoise() > minSignal && peaks[c + 2, j].GetSignalToNoise() > minSignal && peaks[c + 3, j].GetSignalToNoise() > minSignal) count[c, j] = 2;
                                            // Peaks 2 & 3 & 4
                                            else if (peaks[c, j] == null && peaks[c + 1, j].GetSignalToNoise() > minSignal && peaks[c + 2, j].GetSignalToNoise() > minSignal && peaks[c + 3, j].GetSignalToNoise() > minSignal) count[c, j] = 3;
                                            else count[c, j] = -1;
                                        }
                                    }
                                    else if (isotopologues == 6)
                                    {
                                        // 3 peak pattern
                                        if (peakCount[c, j] == 3)
                                        {
                                            if (peaks[c + 3, j] == null)
                                            {
                                                // Peaks 1 & 2
                                                if (peaks[c + 2, j] == null && peaks[c, j].GetSignalToNoise() > minSignal && peaks[c + 1, j].GetSignalToNoise() > minSignal) count[c, j] = 0;
                                                // Peaks 1 & 3
                                                else if (peaks[c + 1, j] == null && peaks[c, j].GetSignalToNoise() > minSignal && peaks[c + 2, j].GetSignalToNoise() > minSignal) count[c, j] = 1;
                                                // Peaks 2 & 3
                                                else if (peaks[c, j] == null && peaks[c + 1, j].GetSignalToNoise() > minSignal && peaks[c + 2, j].GetSignalToNoise() > minSignal) count[c, j] = 3;
                                                else count[c, j] = -1;
                                            }
                                            else if (peaks[c + 2, j] == null)
                                            {
                                                // Peaks 1 & 4
                                                if (peaks[c + 1, j] == null && peaks[c, j].GetSignalToNoise() > minSignal && peaks[c + 3, j].GetSignalToNoise() > minSignal) count[c, j] = 2;
                                                // Peaks 2 & 4
                                                else if (peaks[c, j] == null && peaks[c + 1, j].GetSignalToNoise() > minSignal && peaks[c + 3, j].GetSignalToNoise() > minSignal) count[c, j] = 4;
                                                else count[c, j] = -1;
                                            }
                                            // Peaks 3 & 4
                                            else if (peaks[c + 1, j] == null && peaks[c + 2, j].GetSignalToNoise() > minSignal && peaks[c + 3, j].GetSignalToNoise() > minSignal) count[c, j] = 5;
                                        }
                                        // 4 peak pattern
                                        else if (peakCount[c, j] == 4)
                                        {
                                            // Peaks 1 & 2 & 3
                                            if (peaks[c + 3, j] == null && peaks[c, j].GetSignalToNoise() > minSignal && peaks[c + 1, j].GetSignalToNoise() > minSignal && peaks[c + 2, j].GetSignalToNoise() > minSignal) count[c, j] = 0;
                                            // Peaks 1 & 2 & 4
                                            else if (peaks[c + 2, j] == null && peaks[c, j].GetSignalToNoise() > minSignal && peaks[c + 1, j].GetSignalToNoise() > minSignal && peaks[c + 3, j].GetSignalToNoise() > minSignal) count[c, j] = 1;
                                            // Peaks 1 & 3 & 4
                                            else if (peaks[c + 1, j] == null && peaks[c, j].GetSignalToNoise() > minSignal && peaks[c + 2, j].GetSignalToNoise() > minSignal && peaks[c + 3, j].GetSignalToNoise() > minSignal) count[c, j] = 2;
                                            // Peaks 2 & 3 & 4
                                            else if (peaks[c, j] == null && peaks[c + 1, j].GetSignalToNoise() > minSignal && peaks[c + 2, j].GetSignalToNoise() > minSignal && peaks[c + 3, j].GetSignalToNoise() > minSignal) count[c, j] = 3;
                                            else count[c, j] = -1;
                                        }
                                        // 5 peak pattern
                                        else if (peakCount[c, j] == 5)
                                        {
                                            // Peaks 1 & 2 & 3 & 4
                                            if (peaks[c + 4, j] == null && peaks[c, j].GetSignalToNoise() > minSignal && peaks[c + 1, j].GetSignalToNoise() > minSignal && peaks[c + 2, j].GetSignalToNoise() > minSignal && peaks[c + 3, j].GetSignalToNoise() > minSignal) count[c, j] = 0;
                                            // Peaks 1 & 2 & 3 & 5
                                            else if (peaks[c + 3, j] == null && peaks[c, j].GetSignalToNoise() > minSignal && peaks[c + 1, j].GetSignalToNoise() > minSignal && peaks[c + 2, j].GetSignalToNoise() > minSignal && peaks[c + 4, j].GetSignalToNoise() > minSignal) count[c, j] = 1;
                                            // Peaks 1 & 2 & 4 & 5
                                            else if (peaks[c + 2, j] == null && peaks[c, j].GetSignalToNoise() > minSignal && peaks[c + 1, j].GetSignalToNoise() > minSignal && peaks[c + 3, j].GetSignalToNoise() > minSignal && peaks[c + 4, j].GetSignalToNoise() > minSignal) count[c, j] = 2;
                                            // Peaks 1 & 3 & 4 & 5
                                            else if (peaks[c + 1, j] == null && peaks[c, j].GetSignalToNoise() > minSignal && peaks[c + 2, j].GetSignalToNoise() > minSignal && peaks[c + 3, j].GetSignalToNoise() > minSignal && peaks[c + 4, j].GetSignalToNoise() > minSignal) count[c, j] = 3;
                                            // Peaks 2 & 3 & 4 & 5
                                            else if (peaks[c, j] == null && peaks[c + 1, j].GetSignalToNoise() > minSignal && peaks[c + 2, j].GetSignalToNoise() > minSignal && peaks[c + 3, j].GetSignalToNoise() > minSignal && peaks[c + 4, j].GetSignalToNoise() > minSignal) count[c, j] = 4;
                                            else count[c, j] = -1;
                                        }
                                    }
                                }
                            }
                        }
                    }
                    // Non-NeuCode quantification
                    else
                    {
                        double minSignal = Form1.MINIMUMSN * 1.0;
                        count = new int[1, isotopes];
                        for (int j = 0; j < isotopes; j++)
                        {
                            if (peakCount[0, j] == channels || peakCount[0, j] < parent.peaksNeeded) count[0, j] = 0;
                            else
                            {
                                if (channels == 2)
                                {
                                    // 1 peak pattern
                                    if (peakCount[0, j] == 1)
                                    {
                                        if (peaks[1, j] == null)
                                        {
                                            // Peak 1
                                            if (peaks[0, j].GetSignalToNoise() > minSignal) count[0, j] = 0;
                                            else count[0, j] = -1;
                                        }
                                        else if (peaks[0, j] == null)
                                        {
                                            // Peak 2
                                            if (peaks[1, j].GetSignalToNoise() > minSignal) count[0, j] = 1;
                                            else count[0, j] = -1;
                                        }
                                    }
                                }
                                else if (channels == 3)
                                {
                                    // 1 peak pattern
                                    if (peakCount[0, j] == 1)
                                    {
                                        if (peaks[2, j] == null)
                                        {
                                            // Peak 1
                                            if (peaks[1, j] == null && peaks[0, j].GetSignalToNoise() > minSignal) count[0, j] = 0;
                                            // Peak 2
                                            else if (peaks[0, j] == null && peaks[1, j].GetSignalToNoise() > minSignal) count[0, j] = 1;
                                            else count[0, j] = -1;
                                        }
                                        // Peak 3
                                        else if (peaks[2, j].GetSignalToNoise() > minSignal) count[0, j] = 2;
                                        else count[0, j] = -1;
                                    }
                                    // 2 peak pattern
                                    else if (peakCount[0, j] == 2)
                                    {
                                        // Peaks 1 & 2
                                        if (peaks[2, j] == null && peaks[0, j].GetSignalToNoise() > minSignal && peaks[1, j].GetSignalToNoise() > minSignal) count[0, j] = 0;
                                        // Peaks 1 & 3
                                        else if (peaks[1, j] == null && peaks[0, j].GetSignalToNoise() > minSignal && peaks[2, j].GetSignalToNoise() > minSignal) count[0, j] = 1;
                                        // Peaks 2 & 3
                                        else if (peaks[0, j] == null && peaks[1, j].GetSignalToNoise() > minSignal && peaks[2, j].GetSignalToNoise() > minSignal) count[0, j] = 2;
                                        else count[0, j] = -1;
                                    }
                                }
                            }
                        }
                    }
                }
                return count;
            }
        }
        public int charge; // Charge associated with the parent peptide's best PSM
        public double averageNoise // The average noise level of the pair's detected peaks
        {
            get
            {
                double average = 0;
                int numNoisePeaks = 0;
                double noiseTotal = 0;
                if (peaks != null)
                {
                    for (int i = 0; i < parent.numChannels; i++)
                    {
                        for (int j = 0; j < parent.numIsotopes; j++)
                        {
                            if (peaks[i, j] != null)
                            {
                                numNoisePeaks++;
                                noiseTotal += peaks[i, j].Noise;
                            }
                        }
                    }
                    if (numNoisePeaks > 0)
                    {
                        average = noiseTotal / ((double)numNoisePeaks);
                    }
                }
                return average;
            }
        }
        public double injectionTime { get; set; }
        public double retentionTime { get; set; }
        public int[] isotopeCount
        {
            get
            {
                int channels = parent.numChannels;
                int isotopes = parent.numIsotopes;
                int[] count = new int[channels];
                if (peaks != null)
                {
                    for (int i = 0; i < channels; i++)
                    {
                        int isotopeCount = 0;
                        for (int j = 0; j < isotopes; j++)
                        {
                            if (peaks[i, j] != null)
                            {
                                isotopeCount++;
                            }
                        }
                        count[i] = isotopeCount;
                    }
                }
                return count;
            }
        }

        // Creates a pair associated with a given peptide and scan number within a raw file
        public Pair(PeptideID parent, MSDataFile rawFile, int scanNumber, double injectionTime, double retentionTime)
        {
            this.parent = parent;
            this.scanNumber = scanNumber;
            this.rawFile = rawFile;
            peaks = new ILabeledPeak[parent.numChannels, parent.numIsotopes];
            charge = parent.bestPSMs[rawFile].Charge;
            this.injectionTime = injectionTime;
            this.retentionTime = retentionTime;
        }
        
        public bool checkCoalescence(int Isotope, int ScanNumber, int Channel)
        {
            bool coalescence = true;
            int isotope = Isotope;
            int scanNumber = ScanNumber;
            int channel = Channel;
            ILabeledPeak light;
            ILabeledPeak heavy;
            int channels = parent.numChannels;

            //First, look for pairs in more abundant isotopes, if possible
            if (isotope > 0)
            {
                for (int i = 0; i < channels; i += 2)
                {
                    for (int j = isotope; j >= 0; j--)
                    {
                        light = peaks[i, j];
                        heavy = peaks[i + 1, j];
                        if (light != null && heavy != null && light.X > 0 && heavy.X > 0)
                        {
                            coalescence = false;
                        }
                    }
                }
            }

            if (coalescence)
            {
                //Then, check other isotopes from the same scan, if possible
                if (isotope < Form1.NUMISOTOPES - 1)
                {
                    for (int i = 0; i < channels; i += 2)
                    {
                        for (int j = isotope; j < Form1.NUMISOTOPES; j++)
                        {
                            light = peaks[i, j];
                            heavy = peaks[i + 1, j];
                            if (light != null && heavy != null && light.X > 0 && heavy.X > 0)
                            {
                                coalescence = false;
                            }
                        }
                    }
                }
            }

            if (coalescence)
            {
                //Finally, check prior scans, if necessary, starting with the previous scan
                List<Pair> pairs;
                Pair current;
                bool pairFound = false;
                parent.allPairs.TryGetValue(rawFile, out pairs);
                int k = pairs.Count - 2;
                while (!pairFound && k >= 0)
                {
                    current = pairs[k];
                    for (int i = 0; i < channels; i += 2)
                    {
                        for (int j = 0; j < Form1.NUMISOTOPES; j++)
                        {
                            light = current.peaks[i, j];
                            heavy = current.peaks[i + 1, j];
                            if (light != null && heavy != null && light.X > 0 && heavy.X > 0)
                            {
                                coalescence = false;
                                pairFound = true;
                            }
                        }
                    }
                    k--;
                }
            }
            return coalescence;
        }

        public bool checkIsotopeDistribution(int Cluster)
        {
            bool badDistribution = false;

            int countC12 = countElement("C", parent.sequenceNoMods);
            int countH1 = countElement("H", parent.sequenceNoMods);
            int countN14 = countElement("N", parent.sequenceNoMods);
            int countO16 = countElement("O", parent.sequenceNoMods);
            int countS32 = countElement("S", parent.sequenceNoMods);

            double M1RelAbundance = (countC12 * 1.1) + (countH1 * 0.015) + (countN14 * 0.37);
            double M2RelAbundance = (Math.Pow(countC12 * 1.1, 2)) / 200.0 + (countO16 * 0.2) + (countS32 * 4.21);
            double M2M1RelAbundance = (M2RelAbundance / M1RelAbundance) * 100.0;

            Range<double> M1Range = new Range<double>(M1RelAbundance - (0.5 * M1RelAbundance), M1RelAbundance + (0.5 * M1RelAbundance));
            Range<double> M2Range = new Range<double>(M2RelAbundance - (0.5 * M2RelAbundance), M2RelAbundance + (0.5 * M2RelAbundance));
            Range<double> M2M1Range = new Range<double>(M2M1RelAbundance - (0.5 * M2M1RelAbundance), M2M1RelAbundance + (0.5 * M2M1RelAbundance));

            if (parent.numIsotopes == 3)
            {
                if (parent.numIsotopologues > 1)
                {
                    int channelIndex = Cluster * parent.numIsotopologues;
                    for (int i = channelIndex; i < channelIndex + parent.numIsotopologues; i++)
                    {
                        if (isotopeCount[i] == 0) { }
                        else if (isotopeCount[i] == 1)
                        {
                            if (peaks[i, 2] != null && !complete[Cluster, 2])
                            {
                                peaks[i, 2] = null;
                                badDistribution = true;
                            }
                            else if (peaks[i, 1] != null && !complete[Cluster, 1])
                            {
                                peaks[i, 1] = null;
                                badDistribution = true;
                            }
                        }
                        else if (isotopeCount[i] == 2)
                        {
                            if (peaks[i, 0] == null)
                            {
                                if (!M2M1Range.Contains((peaks[i, 2].GetSignalToNoise() / peaks[i, 1].GetSignalToNoise()) * 100.0))
                                {
                                    peaks[i, 1] = null;
                                    peaks[i, 2] = null;
                                    badDistribution = true;
                                }
                            }
                            else if (peaks[i, 1] == null)
                            {
                                if (!M2Range.Contains((peaks[i, 2].GetSignalToNoise() / peaks[i, 0].GetSignalToNoise()) * 100.0))
                                {
                                    peaks[i, 0] = null;
                                    peaks[i, 2] = null;
                                    badDistribution = true;
                                }
                            }
                            else
                            {
                                if (!M1Range.Contains((peaks[i, 1].GetSignalToNoise() / peaks[i, 0].GetSignalToNoise()) * 100.0))
                                {
                                    peaks[i, 0] = null;
                                    peaks[i, 1] = null;
                                    badDistribution = true;
                                }
                            }
                        }
                        else
                        {
                            if (!M1Range.Contains((peaks[i, 1].GetSignalToNoise() / peaks[i, 0].GetSignalToNoise()) * 100.0))
                            {
                                peaks[i, 1] = null;
                                badDistribution = true;

                                if (!M2Range.Contains((peaks[i, 2].GetSignalToNoise() / peaks[i, 0].GetSignalToNoise()) * 100.0))
                                {
                                    peaks[i, 0] = null;
                                    peaks[i, 2] = null;
                                    badDistribution = true;
                                }
                            }
                            else if (!M2Range.Contains((peaks[i, 2].GetSignalToNoise() / peaks[i, 0].GetSignalToNoise()) * 100.0))
                            {
                                peaks[i, 2] = null;
                                badDistribution = true;
                            }
                        }
                    }

                    channelIndex = Cluster * parent.numIsotopologues;
                    for (int i = channelIndex; i < channelIndex + parent.numIsotopologues; i++)
                    {
                        if (isotopeCount[i] == 0) { }
                        else if (isotopeCount[i] == 1)
                        {
                            if (peaks[i, 2] != null && !complete[Cluster, 2])
                            {
                                peaks[i, 2] = null;
                                badDistribution = true;
                            }
                            else if (peaks[i, 1] != null && !complete[Cluster, 1])
                            {
                                peaks[i, 1] = null;
                                badDistribution = true;
                            }
                        }
                    }
                }
                else
                {
                    int i = Cluster;
                    if (isotopeCount[i] == 0) { }
                    else if (isotopeCount[i] == 1)
                    {
                        if (peaks[i, 2] != null && !complete[0, 2])
                        {
                            peaks[i, 2] = null;
                            badDistribution = true;
                        }
                        else if (peaks[i, 1] != null && !complete[0, 1])
                        {
                            peaks[i, 1] = null;
                            badDistribution = true;
                        }
                    }
                    else if (isotopeCount[i] == 2)
                    {
                        if (peaks[i, 0] == null)
                        {
                            if (!M2M1Range.Contains((peaks[i, 2].GetSignalToNoise() / peaks[i, 1].GetSignalToNoise()) * 100.0))
                            {
                                peaks[i, 1] = null;
                                peaks[i, 2] = null;
                                badDistribution = true;
                            }
                        }
                        else if (peaks[i, 1] == null)
                        {
                            if (!M2Range.Contains((peaks[i, 2].GetSignalToNoise() / peaks[i, 0].GetSignalToNoise()) * 100.0))
                            {
                                peaks[i, 0] = null;
                                peaks[i, 2] = null;
                                badDistribution = true;
                            }
                        }
                        else
                        {
                            if (!M1Range.Contains((peaks[i, 1].GetSignalToNoise() / peaks[i, 0].GetSignalToNoise()) * 100.0))
                            {
                                peaks[i, 0] = null;
                                peaks[i, 1] = null;
                                badDistribution = true;
                            }
                        }
                    }
                    else
                    {
                        if (!M1Range.Contains((peaks[i, 1].GetSignalToNoise() / peaks[i, 0].GetSignalToNoise()) * 100.0))
                        {
                            peaks[i, 1] = null;
                            badDistribution = true;

                            if (!M2Range.Contains((peaks[i, 2].GetSignalToNoise() / peaks[i, 0].GetSignalToNoise()) * 100.0))
                            {
                                peaks[i, 0] = null;
                                peaks[i, 2] = null;
                                badDistribution = true;
                            }
                        }
                        else if (!M2Range.Contains((peaks[i, 2].GetSignalToNoise() / peaks[i, 0].GetSignalToNoise()) * 100.0))
                        {
                            peaks[i, 2] = null;
                            badDistribution = true;
                        }
                    }
                }
            }
            else if (parent.numIsotopes == 2)
            {
                if (parent.numIsotopologues > 1)
                {
                    int channelIndex = Cluster * parent.numIsotopologues;
                    for (int i = channelIndex; i < channelIndex + parent.numIsotopologues; i++)
                    {
                        if (isotopeCount[i] == 0) { }
                        else if (isotopeCount[i] == 1)
                        {
                            if (peaks[i, 1] != null && !complete[Cluster, 1])
                            {
                                peaks[i, 1] = null;
                                badDistribution = true;
                            }
                        }
                        else
                        {
                            if (!M1Range.Contains((peaks[i, 1].GetSignalToNoise() / peaks[i, 0].GetSignalToNoise()) * 100.0))
                            {
                                peaks[i, 0] = null;
                                peaks[i, 1] = null;
                                badDistribution = true;
                            }
                        }
                    }
                }
                else
                {
                    int i = Cluster;
                    if (isotopeCount[i] == 0) { }
                    else if (isotopeCount[i] == 1)
                    {
                        if (peaks[i, 1] != null && !complete[0, 1])
                        {
                            peaks[i, 1] = null;
                            badDistribution = true;
                        }
                    }
                    else
                    {
                        if (!M1Range.Contains((peaks[i, 1].GetSignalToNoise() / peaks[i, 0].GetSignalToNoise()) * 100.0))
                        {
                            peaks[i, 0] = null;
                            peaks[i, 1] = null;
                            badDistribution = true;
                        }
                    }
                }
            }
            return badDistribution;
        }

        public static int countElement(string Element, string Sequence)
        {
            int count = 0;

            Element element = null;
            CSMSL.Chemistry.Element.PeriodicTable.TryGetElement(Element, out element);

            if (element != null)
            {
                for (int i = 0; i < Sequence.Length; i++)
                {
                    AminoAcid aa = AminoAcid.GetResidue(Sequence[i]);
                    count += aa.ChemicalFormula.Count(element);
                }
            }

            return count;
        }

        public double GetMaxChannelIntensity(int channelIndex)
        {
            double maximum = 0;
            for (int j = 0; j < parent.numIsotopes; j++)
            {
                if (peaks[channelIndex, j] != null && peaks[channelIndex, j].GetDenormalizedIntensity(injectionTime) > maximum)
                {
                    maximum = peaks[channelIndex, j].GetDenormalizedIntensity(injectionTime);
                }
            }
            return maximum;
        }

        public bool[,] GetSpacingCheckArray(int Cluster, int Isotope)
        {
            bool[,] spacingCheckArray = null;
            double expSpacing;
            int channelIndex;

            if (parent.numIsotopologues > 1 && parent.numIsotopologueLabels > 0)
            {
                channelIndex = Cluster * parent.numIsotopologues;
                spacingCheckArray = new bool[parent.numIsotopologues, parent.numIsotopologues];

                for (int r = 0; r < parent.numIsotopologues; r++)
                {
                    for (int c = 0; c < parent.numIsotopologues; c++)
                    {
                        if (r == c) spacingCheckArray[r, c] = false;
                        else
                        {
                            if (peaks[r + channelIndex, Isotope] != null && peaks[c + channelIndex, Isotope] != null)
                            {
                                expSpacing = Mass.MassFromMz(peaks[r + channelIndex, Isotope].X, charge) - Mass.MassFromMz(peaks[c + channelIndex, Isotope].X, charge);
                                if (expSpacing < 0) expSpacing = expSpacing * -1.0;
                                spacingCheckArray[r, c] = parent.spacingMassRange[r, c].Contains(expSpacing);
                            }
                            else if (peaks[r + channelIndex, Isotope] != null || peaks[c + channelIndex, Isotope] != null)
                            {
                                spacingCheckArray[r, c] = true;
                            }
                            else
                            {
                                spacingCheckArray[r, c] = false;
                            }
                        }
                    }
                }
            }
            else if (parent.numIsotopologues < 2 && parent.numClusterLabels > 0)
            {
                spacingCheckArray = new bool[parent.numChannels, parent.numChannels];

                for (int r = 0; r < parent.numChannels; r++)
                {
                    for (int c = 0; c < parent.numChannels; c++)
                    {
                        if (r == c) spacingCheckArray[r, c] = false;
                        else
                        {
                            if (peaks[r, Isotope] != null && peaks[c, Isotope] != null)
                            {
                                expSpacing = Mass.MassFromMz(peaks[r, Isotope].X, charge) - Mass.MassFromMz(peaks[c, Isotope].X, charge);
                                if (expSpacing < 0) expSpacing = expSpacing * -1.0;
                                spacingCheckArray[r, c] = parent.spacingMassRange[r, c].Contains(expSpacing);
                            }
                            else if (peaks[r, Isotope] != null || peaks[c, Isotope] != null)
                            {
                                spacingCheckArray[r, c] = true;
                            }
                            else
                            {
                                spacingCheckArray[r, c] = false;
                            }
                        }
                    }
                }
            }
            return spacingCheckArray;
        }
    }
}
