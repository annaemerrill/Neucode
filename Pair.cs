using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Coon;
using CoonThermo.IO;
using Coon.Spectra;

namespace OMNE
{
    class Pair
    {
        public PeptideID parent;
        public int scanNumber;
        public ThermoRawFile rawFile;
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
        }
        //public IsotopePair[] isotopes; //Isotope will be null only if neither channel produces a peak
        public LabeledPeak[,] peaks; //Number of rows equals the number of channels; number of columns equals the number of isotopes considered
        public bool[,] complete
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
                                    if (peaks[i, j] == null || peaks[i, j].MZ == 0)
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
        public double[,] maxIntensity
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
                                    if (peaks[i, j] != null && peaks[i, j].dNL > max[c,j])
                                    {
                                        max[c,j] = peaks[i, j].dNL;
                                    }
                                }
                            }
                        }
                    }
                }
                return max;
            }
        }
        public int[,] peakCount
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

        public Pair(PeptideID parent, ThermoRawFile rawFile, int scanNumber)
        {
            this.parent = parent;
            this.scanNumber = scanNumber;
            this.rawFile = rawFile;
            peaks = new LabeledPeak[parent.numChannels, parent.numIsotopes];
        }
        
        public bool checkCoalescence(int Isotope, int ScanNumber, int Channel)
        {
            bool coalescence = true;
            int isotope = Isotope;
            int scanNumber = ScanNumber;
            int channel = Channel;
            LabeledPeak light;
            LabeledPeak heavy;
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
                        if (light != null && heavy != null && light.MZ > 0 && heavy.MZ > 0)
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
                            if (light != null && heavy != null && light.MZ > 0 && heavy.MZ > 0)
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
                parent.allHILACPairs.TryGetValue(rawFile, out pairs);
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
                            if (light != null && heavy != null && light.MZ > 0 && heavy.MZ > 0)
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
        
        /*public bool checkPeaks(ThermoRawFileScan current)
        {
            bool noiseBandCapped = false;
            int numChannels = Form1.NUMCHANNELS;
            int numIsotopes = Form1.NUMISOTOPES;
            double injectionTime = rawFile.GetInjectionTime(current.ScanNum);
            double noise;
            List<double> noisePeaks = new List<double>();
            PeptideSpectralMatch best;
            int charge;
            LabeledPeak light;
            LabeledPeak heavy;
            Spacing spacing;

            best = parent.bestPSM[rawFile];
            charge = best.Charge; 

            //Calculate the average noise & check for exceeded maximums
            for (int i = 0; i < numChannels; i++)
            {
                for (int j = 0; j < numIsotopes; j++)
                {
                    if (peaks[i, j] != null)
                    {
                        noisePeaks.Add(peaks[i, j].Noise * injectionTime);

                        if (peaks[i, j].Intensity > parent.maxIntensity[i, 0])
                        {
                            parent.maxIntensity[i, 0] = peaks[i, j].Intensity;
                        }
                    }
                }
            }
            noise = noisePeaks.Sum() / (double)noisePeaks.Count;

            //Apply noise to missing channels
            for (int i = 0; i < numChannels; i += 2)
            {
                for (int j = 0; j < numIsotopes; j++)
                {
                    light = peaks[i, j];
                    heavy = peaks[i + 1, j];
                    if (light == null &&  heavy == null)
                    {

                    }
                    else if (light != null && heavy != null)
                    {
                        spacing = new Spacing(light, heavy, charge);
                        parent.peakSpacings.Add(spacing);
                    }
                    else //Check for coalescence
                    {
                        bool coalescence = checkCoalescence(j, current.ScanNum, i);

                        if (coalescence) //Probable coalescence
                        {
                            //Console.WriteLine("possible coalescence detected");
                            if (light == null)
                            {
                                parent.coalescencedNL.Add(heavy.Intensity);
                                parent.coalescenceTIC.Add(heavy.Intensity / (current.TIC * injectionTime));
                            }

                            if (heavy == null)
                            {
                                parent.coalescencedNL.Add(light.Intensity);
                                parent.coalescenceTIC.Add(light.Intensity / (current.TIC * injectionTime));
                            }
                            heavy = null;
                            light = null;
                        }
                        else //No signs of coalescence -- noise-band cap
                        {
                            if (light == null)
                            {
                                light = new LabeledPeak();
                                light.Intensity = noise;
                                light.MZ = 0;
                                noiseBandCapped = true;
                            }
                            if (heavy == null)
                            {
                                heavy = new LabeledPeak();
                                heavy.Intensity = noise;
                                heavy.MZ = 0;
                                noiseBandCapped = true;
                            }
                        }
                    }
                }
            }
            return noiseBandCapped;
        }*/

    }
}
