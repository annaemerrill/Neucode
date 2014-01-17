using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using CSMSL;
using CSMSL.Spectral;
using CSMSL.IO;
using CSMSL.Chemistry;

namespace Coon.NeuQuant
{
    public class IsotopePair
    {
        public Pair Parent;
        public int Isotope;
        public ILabeledPeak[] ChannelPeaks;
        public int TotalPeakCount
        {
            get
            {
                int peakCount = 0;

                for (int i = 0; i < ChannelPeaks.Length; i++)
                {
                    if (ChannelPeaks[i] != null) peakCount++;
                }

                return peakCount;
            }
        }
        public int[] PeakCountByCluster
        {
            get
            {
                int clusters = Parent.Parent.numClusters;
                int isotopologues = Parent.Parent.numIsotopologues;
                if (isotopologues < 2) clusters = 1;
                int[] peakCountByCluster = new int[clusters];
                if (clusters == 1)
                {
                    peakCountByCluster[0] = TotalPeakCount;
                }
                else
                {
                    for (int c = 0; c < clusters; c++)
                    {
                        int channelIndex = c * isotopologues;
                        for (int ch = channelIndex; ch < channelIndex + isotopologues; ch++)
                        {
                            if (ChannelPeaks[ch] != null) peakCountByCluster[c]++;
                        }
                    }
                }
                return peakCountByCluster;
            }
        }
        public bool Complete
        {
            get
            {
                for (int i = 0; i < Parent.Parent.numChannels; i++)
                {
                    if (ChannelPeaks[i] == null) return false;
                    else if (ChannelPeaks[i].X == 0.0) return false;
                }
                return true;
            }
        }
        public bool[] CompleteByCluster
        {
            get
            {
                int clusters = Parent.Parent.numClusters;
                int isotopologues = Parent.Parent.numIsotopologues;
                if (isotopologues < 2) clusters = 1;
                bool[] completeByCluster = new bool[clusters];

                if (clusters == 1)
                {
                    completeByCluster[0] = Complete;
                }
                else
                {
                    for (int c = 0; c < clusters; c++)
                    {
                        bool complete = true;
                        int channelIndex = c * isotopologues;
                        for (int i = channelIndex; i < channelIndex + isotopologues; i++)
                        {
                            if (ChannelPeaks[i] == null) complete = false;
                            else if (ChannelPeaks[i].X == 0.0) complete = false;
                        }
                        completeByCluster[c] = complete;                        
                    }
                }
                return completeByCluster;
            }
        }
        public double MaxIntensity
        {
            get
            {
                double maxIntensity = 0.0;
                for (int c = 0; c < ChannelPeaks.Length; c++)
                {
                    if (ChannelPeaks[c] != null && ChannelPeaks[c].GetDenormalizedIntensity(Parent.InjectionTime) > maxIntensity) maxIntensity = ChannelPeaks[c].GetDenormalizedIntensity(Parent.InjectionTime);
                }
                return maxIntensity;
            }
        }
        public double[] MaxIntensityByCluster
        {
            get
            {
                int clusters = Parent.Parent.numClusters;
                int isotopologues = Parent.Parent.numIsotopologues;
                if (isotopologues < 2) clusters = 1;
                double[] maxByCluster = new double[clusters];

                if (clusters == 1)
                {
                    maxByCluster[0] = MaxIntensity;
                }
                else
                {
                    for (int c = 0; c < clusters; c++)
                    {
                        int channelIndex = c * isotopologues;
                        for (int ch = channelIndex; ch < channelIndex + isotopologues; ch++)
                        {
                            if (ChannelPeaks[ch] != null && ChannelPeaks[ch].GetDenormalizedIntensity(Parent.InjectionTime) > maxByCluster[c]) maxByCluster[c] = ChannelPeaks[ch].GetDenormalizedIntensity(Parent.InjectionTime);
                        }
                    }
                }
                return maxByCluster;
            }

        }
        public double[] MeanNormalizedIntensities
        {
            get
            {
                double[] normalizedIntensities = new double[ChannelPeaks.Length];
                double averageIntensity = 0.0;
                double numPeaks = ChannelPeaks.Length;

                if (TotalPeakCount > 0)
                {
                    for (int c = 0; c < ChannelPeaks.Length; c++)
                    {
                        averageIntensity += ChannelPeaks[c].GetDenormalizedIntensity(Parent.InjectionTime);
                    }
                    averageIntensity = averageIntensity / numPeaks;

                    for (int c = 0; c < ChannelPeaks.Length; c++)
                    {
                        normalizedIntensities[c] = ChannelPeaks[c].GetDenormalizedIntensity(Parent.InjectionTime) / averageIntensity;
                    }
                }
                return normalizedIntensities;
            }
        }
        public int[] PeakPattern
        {
            get
            {
                int channels = Parent.Parent.numChannels;
                int isotopes = Parent.Parent.numIsotopes;
                int clusters = Parent.Parent.numClusters;
                int isotopologues = Parent.Parent.numIsotopologues;
                int[] count = null;
                if (ChannelPeaks != null)
                {
                    // NeuCode quantification
                    if (isotopologues > 1 && clusters > 1)
                    {
                        double minSignal = Form1.MINIMUMSN * 1.0;
                        count = new int[clusters];
                        for (int c = 0; c < clusters; c++)
                        {
                                int channelIndex = c * isotopologues;
                                if (PeakCountByCluster[c] == isotopologues || PeakCountByCluster[c] < Parent.Parent.peaksNeeded) count[c] = 0;
                                else
                                {
                                    if (isotopologues == 2)
                                    {
                                        // 1 peak pattern
                                        if (PeakCountByCluster[c] == 1)
                                        {
                                            if (ChannelPeaks[channelIndex + 1] == null)
                                            {
                                                // Peak 1
                                                if (ChannelPeaks[channelIndex].GetSignalToNoise() > minSignal) count[c] = 0;
                                                else count[c] = -1;
                                            }
                                            else if (ChannelPeaks[channelIndex] == null)
                                            {
                                                // Peak 2
                                                if (ChannelPeaks[channelIndex + 1].GetSignalToNoise() > minSignal) count[c] = 1;
                                                else count[c] = -1;
                                            }
                                        }
                                    }
                                    else if (isotopologues == 3)
                                    {
                                        // 1 peak pattern
                                        if (PeakCountByCluster[c] == 1)
                                        {
                                            if (ChannelPeaks[channelIndex + 2] == null)
                                            {
                                                // Peak 1
                                                if (ChannelPeaks[channelIndex + 1] == null && ChannelPeaks[channelIndex].GetSignalToNoise() > minSignal) count[c] = 0;
                                                // Peak 2
                                                else if (ChannelPeaks[channelIndex] == null && ChannelPeaks[channelIndex + 1].GetSignalToNoise() > minSignal) count[c] = 1;
                                                else count[c] = -1;
                                            }
                                            // Peak 3
                                            else if (ChannelPeaks[channelIndex + 2].GetSignalToNoise() > minSignal) count[c] = 2;
                                            else count[c] = -1;
                                        }
                                        // 2 peak pattern
                                        else if (PeakCountByCluster[c] == 2)
                                        {
                                            // Peaks 1 & 2
                                            if (ChannelPeaks[channelIndex + 2] == null && ChannelPeaks[channelIndex].GetSignalToNoise() > minSignal && ChannelPeaks[channelIndex + 1].GetSignalToNoise() > minSignal) count[c] = 0;
                                            // Peaks 1 & 3
                                            else if (ChannelPeaks[channelIndex + 1] == null && ChannelPeaks[channelIndex].GetSignalToNoise() > minSignal && ChannelPeaks[channelIndex + 2].GetSignalToNoise() > minSignal) count[c] = 1;
                                            // Peaks 2 & 3
                                            else if (ChannelPeaks[channelIndex] == null && ChannelPeaks[channelIndex + 1].GetSignalToNoise() > minSignal && ChannelPeaks[channelIndex + 2].GetSignalToNoise() > minSignal) count[c] = 2;
                                            else count[c] = -1;
                                        }
                                    }
                                    else if (isotopologues == 4)
                                    {
                                        // 2 peak pattern
                                        if (PeakCountByCluster[c] == 2)
                                        {
                                            if (ChannelPeaks[channelIndex + 3] == null)
                                            {
                                                // Peaks 1 & 2
                                                if (ChannelPeaks[channelIndex + 2] == null && ChannelPeaks[channelIndex].GetSignalToNoise() > minSignal && ChannelPeaks[channelIndex + 1].GetSignalToNoise() > minSignal) count[c] = 0;
                                                // Peaks 1 & 3
                                                else if (ChannelPeaks[channelIndex + 1] == null && ChannelPeaks[channelIndex].GetSignalToNoise() > minSignal && ChannelPeaks[channelIndex + 2].GetSignalToNoise() > minSignal) count[c] = 1;
                                                // Peaks 2 & 3
                                                else if (ChannelPeaks[channelIndex] == null && ChannelPeaks[channelIndex + 1].GetSignalToNoise() > minSignal && ChannelPeaks[channelIndex + 2].GetSignalToNoise() > minSignal) count[c] = 3;
                                                else count[c] = -1;
                                            }
                                            else if (ChannelPeaks[channelIndex + 2] == null)
                                            {
                                                // Peaks 1 & 4
                                                if (ChannelPeaks[channelIndex + 1] == null && ChannelPeaks[channelIndex].GetSignalToNoise() > minSignal && ChannelPeaks[channelIndex + 3].GetSignalToNoise() > minSignal) count[c] = 2;
                                                // Peaks 2 & 4
                                                else if (ChannelPeaks[channelIndex] == null && ChannelPeaks[channelIndex + 1].GetSignalToNoise() > minSignal && ChannelPeaks[channelIndex + 3].GetSignalToNoise() > minSignal) count[c] = 4;
                                                else count[c] = -1;
                                            }
                                            // Peaks 3 & 4
                                            else if (ChannelPeaks[channelIndex + 1] == null && ChannelPeaks[channelIndex + 2].GetSignalToNoise() > minSignal && ChannelPeaks[channelIndex + 3].GetSignalToNoise() > minSignal) count[c] = 5;
                                        }
                                        // 3 peak pattern
                                        else if (PeakCountByCluster[c] == 3)
                                        {
                                            // Peaks 1 & 2 & 3
                                            if (ChannelPeaks[channelIndex + 3] == null && ChannelPeaks[channelIndex].GetSignalToNoise() > minSignal && ChannelPeaks[channelIndex + 1].GetSignalToNoise() > minSignal && ChannelPeaks[channelIndex + 2].GetSignalToNoise() > minSignal) count[c] = 0;
                                            // Peaks 1 & 2 & 4
                                            else if (ChannelPeaks[channelIndex + 2] == null && ChannelPeaks[channelIndex].GetSignalToNoise() > minSignal && ChannelPeaks[channelIndex + 1].GetSignalToNoise() > minSignal && ChannelPeaks[channelIndex + 3].GetSignalToNoise() > minSignal) count[c] = 1;
                                            // Peaks 1 & 3 & 4
                                            else if (ChannelPeaks[channelIndex + 1] == null && ChannelPeaks[channelIndex].GetSignalToNoise() > minSignal && ChannelPeaks[channelIndex + 2].GetSignalToNoise() > minSignal && ChannelPeaks[channelIndex + 3].GetSignalToNoise() > minSignal) count[c] = 2;
                                            // Peaks 2 & 3 & 4
                                            else if (ChannelPeaks[channelIndex] == null && ChannelPeaks[channelIndex + 1].GetSignalToNoise() > minSignal && ChannelPeaks[channelIndex + 2].GetSignalToNoise() > minSignal && ChannelPeaks[channelIndex + 3].GetSignalToNoise() > minSignal) count[c] = 3;
                                            else count[c] = -1;
                                        }
                                    }
                                    else if (isotopologues == 6)
                                    {
                                        // 3 peak pattern
                                        if (PeakCountByCluster[c] == 3)
                                        {
                                            // Peaks 1 & 2 & 3
                                            // Peaks 1 & 2 & 4
                                            // Peaks 1 & 2 & 5
                                            // Peaks 1 & 2 & 6
                                            // Peaks 1 & 3 & 4
                                            // Peaks 1 & 3 & 5
                                            // Peaks 1 & 3 & 6
                                            // Peaks 1 & 4 & 5
                                            // Peaks 1 & 4 & 6
                                            // Peaks 1 & 5 & 6
                                            // Peaks 2 & 3 & 4
                                            // Peaks 2 & 3 & 5
                                            // Peaks 2 & 3 & 6
                                            // Peaks 2 & 4 & 5
                                            // Peaks 2 & 4 & 6
                                            // Peaks 2 & 5 & 6
                                            // Peaks 3 & 4 & 5
                                            // Peaks 3 & 4 & 6
                                            // Peaks 3 & 5 & 6
                                            // Peaks 4 & 5 & 6

                                            if (ChannelPeaks[c + 3] == null)
                                            {
                                                // Peaks 1 & 2
                                                if (ChannelPeaks[c + 2] == null && ChannelPeaks[c].GetSignalToNoise() > minSignal && ChannelPeaks[c + 1].GetSignalToNoise() > minSignal) count[c] = 0;
                                                // Peaks 1 & 3
                                                else if (ChannelPeaks[c + 1] == null && ChannelPeaks[c].GetSignalToNoise() > minSignal && ChannelPeaks[c + 2].GetSignalToNoise() > minSignal) count[c] = 1;
                                                // Peaks 2 & 3
                                                else if (ChannelPeaks[c] == null && ChannelPeaks[c + 1].GetSignalToNoise() > minSignal && ChannelPeaks[c + 2].GetSignalToNoise() > minSignal) count[c] = 3;
                                                else count[c] = -1;
                                            }
                                            else if (ChannelPeaks[c + 2] == null)
                                            {
                                                // Peaks 1 & 4
                                                if (ChannelPeaks[c + 1] == null && ChannelPeaks[c].GetSignalToNoise() > minSignal && ChannelPeaks[c + 3].GetSignalToNoise() > minSignal) count[c] = 2;
                                                // Peaks 2 & 4
                                                else if (ChannelPeaks[c] == null && ChannelPeaks[c + 1].GetSignalToNoise() > minSignal && ChannelPeaks[c + 3].GetSignalToNoise() > minSignal) count[c] = 4;
                                                else count[c] = -1;
                                            }
                                            // Peaks 3 & 4
                                            else if (ChannelPeaks[c + 1] == null && ChannelPeaks[c + 2].GetSignalToNoise() > minSignal && ChannelPeaks[c + 3].GetSignalToNoise() > minSignal) count[c] = 5;
                                        }
                                        // 4 peak pattern
                                        else if (PeakCountByCluster[c] == 4)
                                        {
                                            if (ChannelPeaks[c + 5] == null)
                                            {
                                                // Peaks 1 & 2 & 3 & 4
                                                if (ChannelPeaks[c + 4] == null && ChannelPeaks[c].GetSignalToNoise() > minSignal && ChannelPeaks[c + 1].GetSignalToNoise() > minSignal && ChannelPeaks[c + 2].GetSignalToNoise() > minSignal && ChannelPeaks[c + 3].GetSignalToNoise() > minSignal) count[c] = 0;
                                                // Peaks 1 & 2 & 3 & 5
                                                else if (ChannelPeaks[c + 3] == null && ChannelPeaks[c].GetSignalToNoise() > minSignal && ChannelPeaks[c + 1].GetSignalToNoise() > minSignal && ChannelPeaks[c + 2].GetSignalToNoise() > minSignal && ChannelPeaks[c + 4].GetSignalToNoise() > minSignal) count[c] = 1;
                                                // Peaks 1 & 2 & 3 & 6
                                                else if (ChannelPeaks[c + 2] == null && ChannelPeaks[c].GetSignalToNoise() > minSignal && ChannelPeaks[c + 1].GetSignalToNoise() > minSignal && ChannelPeaks[c + 3].GetSignalToNoise() > minSignal && ChannelPeaks[c + 4].GetSignalToNoise() > minSignal) count[c] = 2;
                                                // Peaks 1 & 2 & 4 & 5
                                                else if (ChannelPeaks[c + 2] == null && ChannelPeaks[c].GetSignalToNoise() > minSignal && ChannelPeaks[c + 1].GetSignalToNoise() > minSignal && ChannelPeaks[c + 3].GetSignalToNoise() > minSignal && ChannelPeaks[c + 4].GetSignalToNoise() > minSignal) count[c] = 2;
                                            }
                                            else if (ChannelPeaks[c + 4] == null)
                                            {

                                            }
                                            else if (ChannelPeaks[c + 3] == null)
                                            {

                                            }
                                            else if (ChannelPeaks[c + 2] == null)
                                            {

                                            }
                                            else if (ChannelPeaks[c + 1] == null)
                                            {

                                            }
                                            else count[c] = -1;




                                            // Peaks 1 & 2 & 4 & 6
                                            // Peaks 1 & 2 & 5 & 6
                                            // Peaks 1 & 3 & 4 & 5
                                            // Peaks 1 & 3 & 4 & 6
                                            // Peaks 1 & 3 & 5 & 6
                                            // Peaks 1 & 4 & 5 & 6
                                            // Peaks 2 & 3 & 4 & 5
                                            // Peaks 2 & 3 & 4 & 6
                                            // Peaks 2 & 3 & 5 & 6
                                            // Peaks 2 & 4 & 5 & 6
                                            // Peaks 3 & 4 & 5 & 6
                                        }
                                        // 5 peak pattern
                                        else if (PeakCountByCluster[c] == 5)
                                        {
                                            // Peaks 1 & 2 & 3 & 4 & 5
                                            if (ChannelPeaks[c + 5] == null && ChannelPeaks[c].GetSignalToNoise() > minSignal && ChannelPeaks[c + 1].GetSignalToNoise() > minSignal && ChannelPeaks[c + 2].GetSignalToNoise() > minSignal && ChannelPeaks[c + 3].GetSignalToNoise() > minSignal && ChannelPeaks[c + 4].GetSignalToNoise() > minSignal) count[c] = 0;
                                            // Peaks 1 & 2 & 3 & 4 & 6
                                            else if (ChannelPeaks[c + 4] == null && ChannelPeaks[c].GetSignalToNoise() > minSignal && ChannelPeaks[c + 1].GetSignalToNoise() > minSignal && ChannelPeaks[c + 2].GetSignalToNoise() > minSignal && ChannelPeaks[c + 3].GetSignalToNoise() > minSignal && ChannelPeaks[c + 5].GetSignalToNoise() > minSignal) count[c] = 0;
                                            // Peaks 1 & 2 & 3 & 5 & 6
                                            else if (ChannelPeaks[c + 3] == null && ChannelPeaks[c].GetSignalToNoise() > minSignal && ChannelPeaks[c + 1].GetSignalToNoise() > minSignal && ChannelPeaks[c + 2].GetSignalToNoise() > minSignal && ChannelPeaks[c + 4].GetSignalToNoise() > minSignal && ChannelPeaks[c + 5].GetSignalToNoise() > minSignal) count[c] = 0;
                                            // Peaks 1 & 2 & 4 & 6 & 6
                                            else if (ChannelPeaks[c + 2] == null && ChannelPeaks[c].GetSignalToNoise() > minSignal && ChannelPeaks[c + 1].GetSignalToNoise() > minSignal && ChannelPeaks[c + 3].GetSignalToNoise() > minSignal && ChannelPeaks[c + 4].GetSignalToNoise() > minSignal && ChannelPeaks[c + 5].GetSignalToNoise() > minSignal) count[c] = 0;
                                            // Peaks 1 & 3 & 4 & 5 & 6
                                            else if (ChannelPeaks[c + 1] == null && ChannelPeaks[c].GetSignalToNoise() > minSignal && ChannelPeaks[c + 2].GetSignalToNoise() > minSignal && ChannelPeaks[c + 3].GetSignalToNoise() > minSignal && ChannelPeaks[c + 4].GetSignalToNoise() > minSignal && ChannelPeaks[c + 5].GetSignalToNoise() > minSignal) count[c] = 0;
                                            // Peaks 2 & 3 & 4 & 5 & 6
                                            else if (ChannelPeaks[c] == null && ChannelPeaks[c + 1].GetSignalToNoise() > minSignal && ChannelPeaks[c + 2].GetSignalToNoise() > minSignal && ChannelPeaks[c + 3].GetSignalToNoise() > minSignal && ChannelPeaks[c + 4].GetSignalToNoise() > minSignal && ChannelPeaks[c + 5].GetSignalToNoise() > minSignal) count[c] = 0;
                                            else count[c] = -1;
                                        }
                                    }
                                }
                            }
                    }
                    // Non-NeuCode quantification
                    else
                    {
                        double minSignal = Form1.MINIMUMSN * 1.0;
                        count = new int[1];
                        for (int j = 0; j < isotopes; j++)
                        {
                            if (PeakCountByCluster[0] == channels || PeakCountByCluster[0] < Parent.Parent.peaksNeeded) count[0] = 0;
                            else
                            {
                                if (channels == 2)
                                {
                                    // 1 peak pattern
                                    if (PeakCountByCluster[0] == 1)
                                    {
                                        if (ChannelPeaks[1] == null)
                                        {
                                            // Peak 1
                                            if (ChannelPeaks[0].GetSignalToNoise() > minSignal) count[0] = 0;
                                            else count[0] = -1;
                                        }
                                        else if (ChannelPeaks[0] == null)
                                        {
                                            // Peak 2
                                            if (ChannelPeaks[1].GetSignalToNoise() > minSignal) count[0] = 1;
                                            else count[0] = -1;
                                        }
                                    }
                                }
                                else if (channels == 3)
                                {
                                    // 1 peak pattern
                                    if (PeakCountByCluster[0] == 1)
                                    {
                                        if (ChannelPeaks[2] == null)
                                        {
                                            // Peak 1
                                            if (ChannelPeaks[1] == null && ChannelPeaks[0].GetSignalToNoise() > minSignal) count[0] = 0;
                                            // Peak 2
                                            else if (ChannelPeaks[0] == null && ChannelPeaks[1].GetSignalToNoise() > minSignal) count[0] = 1;
                                            else count[0] = -1;
                                        }
                                        // Peak 3
                                        else if (ChannelPeaks[2].GetSignalToNoise() > minSignal) count[0] = 2;
                                        else count[0] = -1;
                                    }
                                    // 2 peak pattern
                                    else if (PeakCountByCluster[0] == 2)
                                    {
                                        // Peaks 1 & 2
                                        if (ChannelPeaks[2] == null && ChannelPeaks[0].GetSignalToNoise() > minSignal && ChannelPeaks[1].GetSignalToNoise() > minSignal) count[0] = 0;
                                        // Peaks 1 & 3
                                        else if (ChannelPeaks[1] == null && ChannelPeaks[0].GetSignalToNoise() > minSignal && ChannelPeaks[2].GetSignalToNoise() > minSignal) count[0] = 1;
                                        // Peaks 2 & 3
                                        else if (ChannelPeaks[0] == null && ChannelPeaks[1].GetSignalToNoise() > minSignal && ChannelPeaks[2].GetSignalToNoise() > minSignal) count[0] = 2;
                                        else count[0] = -1;
                                    }
                                }
                            }
                        }
                    }
                }
                return count;
            }
        }

        // 2 isotopologue pairs
        public IsotopePair(Pair parent, int isotope)
        {
            Parent = parent;
            Isotope = isotope;
            ChannelPeaks = new ILabeledPeak[Parent.Parent.numChannels];
        }

        public bool[,] GetSpacingCheckArray(int Cluster = 0)
        {
            bool[,] spacingCheckArray = null;
            int isotopologues = Parent.Parent.numIsotopologues;
            int channels = Parent.Parent.numChannels;
            double expSpacing;
            int channelIndex;

            if (isotopologues > 1 && Parent.Parent.isotopologueLabel)
            {
                channelIndex = Cluster * isotopologues;
                spacingCheckArray = new bool[isotopologues, isotopologues];

                for (int r = 0; r < isotopologues; r++)
                {
                    for (int c = 0; c < isotopologues; c++)
                    {
                        if (r == c) spacingCheckArray[r, c] = false;
                        else
                        {
                            if (ChannelPeaks[r + channelIndex] != null && ChannelPeaks[c + channelIndex] != null)
                            {
                                expSpacing = Mass.MassFromMz(ChannelPeaks[r + channelIndex].X, Parent.Charge) - Mass.MassFromMz(ChannelPeaks[c + channelIndex].X, Parent.Charge);
                                if (expSpacing < 0) expSpacing = expSpacing * -1.0;
                                spacingCheckArray[r, c] = Parent.Parent.spacingMassRange[r, c].Contains(expSpacing);
                            }
                            else if (ChannelPeaks[r + channelIndex] != null || ChannelPeaks[c + channelIndex] != null)
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
            else if (isotopologues < 2 && Parent.Parent.clusterLabel)
            {
                spacingCheckArray = new bool[channels, channels];

                for (int r = 0; r < channels; r++)
                {
                    for (int c = 0; c < channels; c++)
                    {
                        if (r == c) spacingCheckArray[r, c] = false;
                        else
                        {
                            if (ChannelPeaks[r] != null && ChannelPeaks[c] != null)
                            {
                                expSpacing = Mass.MassFromMz(ChannelPeaks[r].X, Parent.Charge) - Mass.MassFromMz(ChannelPeaks[c].X, Parent.Charge);
                                if (expSpacing < 0) expSpacing = expSpacing * -1.0;
                                spacingCheckArray[r, c] = Parent.Parent.spacingMassRange[r, c].Contains(expSpacing);
                            }
                            else if (ChannelPeaks[r] != null || ChannelPeaks[c] != null)
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
