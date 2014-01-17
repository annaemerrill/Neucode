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


namespace Coon.NeuQuant
{
    public class Pair
    {
        public PeptideID Parent; // The peptide the pair belongs to
        public int ScanNumber; // The scan number the pair is associated with
        public MSDataFile RawFile; // The raw file the pair is associated with
        public int TotalPeakCount
        {
            get
            {
                int count = 0;
                int channels = Parent.numChannels;
                int isotopes = Parent.numIsotopes;
                int clusters = Parent.numClusters;
                if (PeakCount != null)
                {
                    if (Parent.numIsotopologues > 1)
                    {
                        for (int c = 0; c < clusters; c++)
                        {
                            for (int i = 0; i < isotopes; i++)
                            {
                                count += PeakCount[c, i];
                            }
                        }
                    }
                    else
                    {
                        for (int i = 0; i < isotopes; i++)
                        {
                            count += PeakCount[0, i];
                        }
                    }
                }
                return count;
            }
        } // The total number of peaks associated with this pair
        //public ILabeledPeak[,] Peaks; // Keeps track of peaks for quantitation (rows = # channels; columns = # isotopes)
        //public bool[,] Complete // Keeps track of whether or not a cluster (i.e., set of isotopologues) of peaks is complete (rows = # clusters; columns = # isotopes)
        //{
        //    get
        //    {
        //        int channels = parent.numChannels;
        //        int isotopes = parent.numIsotopes;
        //        int clusters = parent.numClusters;
        //        int isotopologues = parent.numIsotopologues;
        //        int channelIndex;
        //        bool[,] complete;
        //        if (isotopologues > 1)
        //        {
        //            complete = new bool[clusters, isotopes];
        //            if (peaks != null)
        //            {
        //                for (int c = 0; c < clusters; c++)
        //                {
        //                    for (int j = 0; j < isotopes; j++)
        //                    {
        //                        complete[c, j] = true;
        //                        channelIndex = c * isotopologues;
        //                        for (int i = channelIndex; i < channelIndex + isotopologues; i++)
        //                        {
        //                            // A cluster is incomplete if any of its peaks are null or noise-band capped
        //                            if (peaks[i, j] == null || peaks[i, j].X == 0)
        //                            {
        //                                complete[c, j] = false;
        //                            }
        //                        }
        //                    }
        //                }
        //            }
        //            return complete;
        //        }
        //        else
        //        {
        //            complete = new bool[1, isotopes];
        //            if (peaks != null)
        //            {
        //                for (int j = 0; j < isotopes; j++)
        //                {
        //                    complete[0,j] = true;
        //                    for (int i = 0; i < channels; i++)
        //                    {
        //                        if (peaks[i, j] == null || peaks[i, j].X == 0)
        //                        {
        //                            complete[0, j] = false;
        //                        }
        //                    }
        //                }
        //            }
        //            return complete;
        //        }
        //    }
        //}
        public double[] MaxClusterIntensity // Keeps track of each cluster's maximum intensity (rows = # clusters; columns = # isotopes)
        {
            get
            {
                int channels = Parent.numChannels;
                int clusters = Parent.numClusters;
                int isotopes = Parent.numIsotopes;
                int isotopologues = Parent.numIsotopologues;
                double[] max = new double[clusters];
                {
                    if (IsotopePairs.Count > 0)
                    {
                        foreach (IsotopePair isotopePair in IsotopePairs)
                        {
                            for (int c = 0; c < clusters; c++)
                            {
                                if (isotopePair.MaxIntensityByCluster[c] > max[c]) max[c] = isotopePair.MaxIntensityByCluster[c];
                            }
                        }
                    }
                }
                return max;
            }
        }
        public int[,] PeakCount // Keeps track of the number of peaks in each cluster (rows = # clusters; columns = # isotopes)
        {
            get
            {
                int channels = Parent.numChannels;
                int isotopes = Parent.numIsotopes;
                int clusters = Parent.numClusters;
                int isotopologues = Parent.numIsotopologues;
                if (isotopologues < 2) clusters = 1;
                int[,] count = new int[clusters, isotopes];
                if (IsotopePairs.Count > 0)
                {
                    foreach (IsotopePair isotopePair in IsotopePairs)
                    {
                        for (int c = 0; c < clusters; c++)
                        {
                            count[c, isotopePair.Isotope] = isotopePair.PeakCountByCluster[c];
                        }
                    }
                }
                return count;
            }
        }        
        public int Charge; // Charge associated with the parent peptide's best PSM
        public double AverageNoise // The average noise level of the pair's detected peaks
        {
            get
            {
                double average = 0;
                int numNoisePeaks = 0;
                double noiseTotal = 0;
                if (IsotopePairs.Count > 0)
                {
                    foreach (IsotopePair isotopePair in IsotopePairs)
                    {
                        for (int i = 0; i < Parent.numChannels; i++)
                        {
                            if (isotopePair.ChannelPeaks[i] != null)
                            {
                                numNoisePeaks++;
                                noiseTotal += isotopePair.ChannelPeaks[i].Noise;
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
        public double InjectionTime { get; set; }
        public double RetentionTime { get; set; }
        public int[] IsotopeCount
        {
            get
            {
                int channels = Parent.numChannels;
                int isotopes = Parent.numIsotopes;
                int[] count = new int[channels];
                if (IsotopePairs.Count > 0)
                {        
                    for (int i = 0; i < channels; i++)
                    {
                        int isotopeCount = 0;
                        foreach (IsotopePair isotopePair in IsotopePairs)
                        {
                            if (isotopePair.ChannelPeaks[i] != null)
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
        public int[] PeaksPerIsotope
        {
            get
            {
                int channels = Parent.numChannels;
                int isotopes = Parent.numIsotopes;
                int[] peaksPerIsotope = new int[isotopes];

                foreach (IsotopePair isotopePair in IsotopePairs)
                {
                    if (isotopePair != null) peaksPerIsotope[isotopePair.Isotope] = isotopePair.TotalPeakCount;
                }

                return peaksPerIsotope;
            }
        }
        public List<IsotopePair> IsotopePairs;

        // Creates a pair associated with a given peptide and scan number within a raw file
        public Pair(PeptideID parent, MSDataFile rawFile, int scanNumber, int charge, double injectionTime, double retentionTime)
        {
            Parent = parent;
            ScanNumber = scanNumber;
            RawFile = rawFile;
            Charge = charge;
            InjectionTime = injectionTime;
            RetentionTime = retentionTime;
            IsotopePairs = new List<IsotopePair>();
        }
        
        //public bool checkCoalescence(int Isotope, int ScanNumber, int Channel)
        //{
        //    bool coalescence = true;
        //    int isotope = Isotope;
        //    int scanNumber = ScanNumber;
        //    int channel = Channel;
        //    ILabeledPeak light;
        //    ILabeledPeak heavy;
        //    int channels = parent.numChannels;

        //    //First, look for pairs in more abundant isotopes, if possible
        //    if (isotope > 0)
        //    {
        //        for (int i = 0; i < channels; i += 2)
        //        {
        //            for (int j = isotope; j >= 0; j--)
        //            {
        //                light = peaks[i, j];
        //                heavy = peaks[i + 1, j];
        //                if (light != null && heavy != null && light.X > 0 && heavy.X > 0)
        //                {
        //                    coalescence = false;
        //                }
        //            }
        //        }
        //    }

        //    if (coalescence)
        //    {
        //        //Then, check other isotopes from the same scan, if possible
        //        if (isotope < Form1.NUMISOTOPES - 1)
        //        {
        //            for (int i = 0; i < channels; i += 2)
        //            {
        //                for (int j = isotope; j < Form1.NUMISOTOPES; j++)
        //                {
        //                    light = peaks[i, j];
        //                    heavy = peaks[i + 1, j];
        //                    if (light != null && heavy != null && light.X > 0 && heavy.X > 0)
        //                    {
        //                        coalescence = false;
        //                    }
        //                }
        //            }
        //        }
        //    }

        //    if (coalescence)
        //    {
        //        //Finally, check prior scans, if necessary, starting with the previous scan
        //        List<Pair> pairs;
        //        Pair current;
        //        bool pairFound = false;
        //        parent.allPairs.TryGetValue(rawFile, out pairs);
        //        int k = pairs.Count - 2;
        //        while (!pairFound && k >= 0)
        //        {
        //            current = pairs[k];
        //            for (int i = 0; i < channels; i += 2)
        //            {
        //                for (int j = 0; j < Form1.NUMISOTOPES; j++)
        //                {
        //                    light = current.peaks[i, j];
        //                    heavy = current.peaks[i + 1, j];
        //                    if (light != null && heavy != null && light.X > 0 && heavy.X > 0)
        //                    {
        //                        coalescence = false;
        //                        pairFound = true;
        //                    }
        //                }
        //            }
        //            k--;
        //        }
        //    }
        //    return coalescence;
        //}

        //public bool checkIsotopeDistribution(int Cluster)
        //{
        //    bool badDistribution = false;

        //    int countC12 = countElement("C", parent.sequenceNoMods);
        //    int countH1 = countElement("H", parent.sequenceNoMods);
        //    int countN14 = countElement("N", parent.sequenceNoMods);
        //    int countO16 = countElement("O", parent.sequenceNoMods);
        //    int countS32 = countElement("S", parent.sequenceNoMods);

        //    double M1RelAbundance = (countC12 * 1.1) + (countH1 * 0.015) + (countN14 * 0.37);
        //    double M2RelAbundance = (Math.Pow(countC12 * 1.1, 2)) / 200.0 + (countO16 * 0.2) + (countS32 * 4.21);
        //    double M2M1RelAbundance = (M2RelAbundance / M1RelAbundance) * 100.0;

        //    Range<double> M1Range = new Range<double>(M1RelAbundance - (0.5 * M1RelAbundance), M1RelAbundance + (0.5 * M1RelAbundance));
        //    Range<double> M2Range = new Range<double>(M2RelAbundance - (0.5 * M2RelAbundance), M2RelAbundance + (0.5 * M2RelAbundance));
        //    Range<double> M2M1Range = new Range<double>(M2M1RelAbundance - (0.5 * M2M1RelAbundance), M2M1RelAbundance + (0.5 * M2M1RelAbundance));

        //    if (parent.numIsotopes == 3)
        //    {
        //        if (parent.numIsotopologues > 1)
        //        {
        //            int channelIndex = Cluster * parent.numIsotopologues;
        //            for (int i = channelIndex; i < channelIndex + parent.numIsotopologues; i++)
        //            {
        //                if (isotopeCount[i] == 0) { }
        //                else if (isotopeCount[i] == 1)
        //                {
        //                    if (peaks[i, 2] != null && !complete[Cluster, 2])
        //                    {
        //                        peaks[i, 2] = null;
        //                        badDistribution = true;
        //                    }
        //                    else if (peaks[i, 1] != null && !complete[Cluster, 1])
        //                    {
        //                        peaks[i, 1] = null;
        //                        badDistribution = true;
        //                    }
        //                }
        //                else if (isotopeCount[i] == 2)
        //                {
        //                    if (peaks[i, 0] == null)
        //                    {
        //                        if (!M2M1Range.Contains((peaks[i, 2].GetSignalToNoise() / peaks[i, 1].GetSignalToNoise()) * 100.0))
        //                        {
        //                            peaks[i, 1] = null;
        //                            peaks[i, 2] = null;
        //                            badDistribution = true;
        //                        }
        //                    }
        //                    else if (peaks[i, 1] == null)
        //                    {
        //                        if (!M2Range.Contains((peaks[i, 2].GetSignalToNoise() / peaks[i, 0].GetSignalToNoise()) * 100.0))
        //                        {
        //                            peaks[i, 0] = null;
        //                            peaks[i, 2] = null;
        //                            badDistribution = true;
        //                        }
        //                    }
        //                    else
        //                    {
        //                        if (!M1Range.Contains((peaks[i, 1].GetSignalToNoise() / peaks[i, 0].GetSignalToNoise()) * 100.0))
        //                        {
        //                            peaks[i, 0] = null;
        //                            peaks[i, 1] = null;
        //                            badDistribution = true;
        //                        }
        //                    }
        //                }
        //                else
        //                {
        //                    if (!M1Range.Contains((peaks[i, 1].GetSignalToNoise() / peaks[i, 0].GetSignalToNoise()) * 100.0))
        //                    {
        //                        peaks[i, 1] = null;
        //                        badDistribution = true;

        //                        if (!M2Range.Contains((peaks[i, 2].GetSignalToNoise() / peaks[i, 0].GetSignalToNoise()) * 100.0))
        //                        {
        //                            peaks[i, 0] = null;
        //                            peaks[i, 2] = null;
        //                            badDistribution = true;
        //                        }
        //                    }
        //                    else if (!M2Range.Contains((peaks[i, 2].GetSignalToNoise() / peaks[i, 0].GetSignalToNoise()) * 100.0))
        //                    {
        //                        peaks[i, 2] = null;
        //                        badDistribution = true;
        //                    }
        //                }
        //            }

        //            channelIndex = Cluster * parent.numIsotopologues;
        //            for (int i = channelIndex; i < channelIndex + parent.numIsotopologues; i++)
        //            {
        //                if (isotopeCount[i] == 0) { }
        //                else if (isotopeCount[i] == 1)
        //                {
        //                    if (peaks[i, 2] != null && !complete[Cluster, 2])
        //                    {
        //                        peaks[i, 2] = null;
        //                        badDistribution = true;
        //                    }
        //                    else if (peaks[i, 1] != null && !complete[Cluster, 1])
        //                    {
        //                        peaks[i, 1] = null;
        //                        badDistribution = true;
        //                    }
        //                }
        //            }
        //        }
        //        else
        //        {
        //            int i = Cluster;
        //            if (isotopeCount[i] == 0) { }
        //            else if (isotopeCount[i] == 1)
        //            {
        //                if (peaks[i, 2] != null && !complete[0, 2])
        //                {
        //                    peaks[i, 2] = null;
        //                    badDistribution = true;
        //                }
        //                else if (peaks[i, 1] != null && !complete[0, 1])
        //                {
        //                    peaks[i, 1] = null;
        //                    badDistribution = true;
        //                }
        //            }
        //            else if (isotopeCount[i] == 2)
        //            {
        //                if (peaks[i, 0] == null)
        //                {
        //                    if (!M2M1Range.Contains((peaks[i, 2].GetSignalToNoise() / peaks[i, 1].GetSignalToNoise()) * 100.0))
        //                    {
        //                        peaks[i, 1] = null;
        //                        peaks[i, 2] = null;
        //                        badDistribution = true;
        //                    }
        //                }
        //                else if (peaks[i, 1] == null)
        //                {
        //                    if (!M2Range.Contains((peaks[i, 2].GetSignalToNoise() / peaks[i, 0].GetSignalToNoise()) * 100.0))
        //                    {
        //                        peaks[i, 0] = null;
        //                        peaks[i, 2] = null;
        //                        badDistribution = true;
        //                    }
        //                }
        //                else
        //                {
        //                    if (!M1Range.Contains((peaks[i, 1].GetSignalToNoise() / peaks[i, 0].GetSignalToNoise()) * 100.0))
        //                    {
        //                        peaks[i, 0] = null;
        //                        peaks[i, 1] = null;
        //                        badDistribution = true;
        //                    }
        //                }
        //            }
        //            else
        //            {
        //                if (!M1Range.Contains((peaks[i, 1].GetSignalToNoise() / peaks[i, 0].GetSignalToNoise()) * 100.0))
        //                {
        //                    peaks[i, 1] = null;
        //                    badDistribution = true;

        //                    if (!M2Range.Contains((peaks[i, 2].GetSignalToNoise() / peaks[i, 0].GetSignalToNoise()) * 100.0))
        //                    {
        //                        peaks[i, 0] = null;
        //                        peaks[i, 2] = null;
        //                        badDistribution = true;
        //                    }
        //                }
        //                else if (!M2Range.Contains((peaks[i, 2].GetSignalToNoise() / peaks[i, 0].GetSignalToNoise()) * 100.0))
        //                {
        //                    peaks[i, 2] = null;
        //                    badDistribution = true;
        //                }
        //            }
        //        }
        //    }
        //    else if (parent.numIsotopes == 2)
        //    {
        //        if (parent.numIsotopologues > 1)
        //        {
        //            int channelIndex = Cluster * parent.numIsotopologues;
        //            for (int i = channelIndex; i < channelIndex + parent.numIsotopologues; i++)
        //            {
        //                if (isotopeCount[i] == 0) { }
        //                else if (isotopeCount[i] == 1)
        //                {
        //                    if (peaks[i, 1] != null && !complete[Cluster, 1])
        //                    {
        //                        peaks[i, 1] = null;
        //                        badDistribution = true;
        //                    }
        //                }
        //                else
        //                {
        //                    if (!M1Range.Contains((peaks[i, 1].GetSignalToNoise() / peaks[i, 0].GetSignalToNoise()) * 100.0))
        //                    {
        //                        peaks[i, 0] = null;
        //                        peaks[i, 1] = null;
        //                        badDistribution = true;
        //                    }
        //                }
        //            }
        //        }
        //        else
        //        {
        //            int i = Cluster;
        //            if (isotopeCount[i] == 0) { }
        //            else if (isotopeCount[i] == 1)
        //            {
        //                if (peaks[i, 1] != null && !complete[0, 1])
        //                {
        //                    peaks[i, 1] = null;
        //                    badDistribution = true;
        //                }
        //            }
        //            else
        //            {
        //                if (!M1Range.Contains((peaks[i, 1].GetSignalToNoise() / peaks[i, 0].GetSignalToNoise()) * 100.0))
        //                {
        //                    peaks[i, 0] = null;
        //                    peaks[i, 1] = null;
        //                    badDistribution = true;
        //                }
        //            }
        //        }
        //    }
        //    return badDistribution;
        //}

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
            foreach (IsotopePair isotopePair in IsotopePairs)
            {
                if (isotopePair.ChannelPeaks[channelIndex] != null && isotopePair.ChannelPeaks[channelIndex].X > 0.0 && isotopePair.ChannelPeaks[channelIndex].GetDenormalizedIntensity(InjectionTime) > maximum)
                {
                    maximum = isotopePair.ChannelPeaks[channelIndex].GetDenormalizedIntensity(InjectionTime);
                }
            }
            return maximum;
        }
    }
}
