using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using CSMSL;
using CSMSL.IO;
using CSMSL.Spectral;


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
                    for (int c = 0; c < clusters; c++)
                    {
                        for (int i = 0; i < isotopes; i++)
                        {
                            count += peakCount[c,i];
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
                bool[,] complete = new bool[clusters, isotopes];
                {
                    if (peaks != null)
                    {
                        for (int c = 0; c < clusters; c++)
                        {
                            for (int j = 0; j < isotopes; j++)
                            {
                                complete[c,j] = true;
                                channelIndex = c * isotopologues;
                                for (int i = channelIndex; i < channelIndex + isotopologues; i++)
                                {
                                    // A cluster is incomplete if any of its peaks are null or noise-band capped
                                    if (peaks[i, j] == null || peaks[i, j].X == 0)
                                    {
                                        complete[c,j] = false;
                                    }
                                }
                            }
                        }
                        return complete;
                    }
                }
                return complete;
            }
        }
        public double[,] maxIntensity // Keeps track of each cluster's maximum intensity (rows = # clusters; columns = # isotopes)
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
                        double injectionTime = rawFile.GetInjectionTime(scanNumber);
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
                int[,] count = new int[clusters, isotopes];
                if (peaks != null)
                {
                    for (int j = 0; j < isotopes; j++)
                    {
                        for (int c = 0; c < clusters; c++)
                        {
                            int channelIndex = c * isotopologues;
                            for (int i = channelIndex; i < channelIndex + isotopologues; i++)
                            {
                                if (peaks[i, j] != null)
                                {
                                    count[c,j]++;
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

        // Creates a pair associated with a given peptide and scan number within a raw file
        public Pair(PeptideID parent, MSDataFile rawFile, int scanNumber)
        {
            this.parent = parent;
            this.scanNumber = scanNumber;
            this.rawFile = rawFile;
            peaks = new ILabeledPeak[parent.numChannels, parent.numIsotopes];
            charge = parent.PSMs[rawFile][0].Charge;
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
    }
}
