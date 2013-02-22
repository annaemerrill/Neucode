using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Windows.Forms;
using System.IO;
using LumenWorks.Framework.IO.Csv;
using CSMSL;
using CSMSL.IO;
using CSMSL.IO.Thermo;
using CSMSL.Spectral;
using CSMSL.Proteomics;
using CSMSL.Chemistry;

namespace Coon.NeuQuant
{
    public partial class Form1 : Form
    {
        public static bool NEUCODE_DUPLEX_LYS;
        public static bool NEUCODE_DUPLEX_LYS1;
        public static bool NEUCODE_DUPLEX_CARBAMYL;
        public static bool NEUCODE_SIXPLEX_MTRAQ;
        public static bool SILAC_DUPLEX_LYS;
        public static bool NEUCODE_SIXPLEX_ARG;
        public static bool NEUCODE_12PLEX;
        public static bool NEUCODE;
        public static bool TRADITIONAL;
        public static double RTWINDOW;
        public static int NUMCHANNELS;
        public static int NUMISOTOPES;
        public static int NUMISOTOPOLOGUES;
        public static int NUMCLUSTERS;
        public static MSDataFile RAWFILE;
        public static double SYSTEMATICERROR;
        public static bool FIRSTSEARCHDONE;
        public static bool CHECKPAIRSPACING;
        public static bool NOISEBANDCAP;
        public static bool PEAKCOALESCENCE;
        public static bool QUANTFILTER;
        public static bool CORRECTPROLINE;
        public static double MINIMUMSN;
        public static double MAXIMUMDNL;
        public static Dictionary<string, MSDataFile> RAWFILES;
        public static bool NEUCODE_ARG_PROLINECONVERSION;
        public static bool MULTIINJECT;
        public static bool CALCIUM;
        public static bool NEUCODE_4PLEX_HEAVY;
        public static bool NEUCODE_4PLEX_MEDIUM;
        public static bool NEUCODE_4PLEX_LIGHT;
        public static Dictionary<int, List<double>> PEAKSNPAIRS;
        public static int QUANTCOUNT;
        
        public Form1()
        {

            InitializeComponent();
            rawFileBox.Text = @"C:\Users\Anna\Desktop\NeuQuant CSMSL Test";
            csvInputBox.Text = @"C:\Users\Anna\Desktop\NeuQuant CSMSL Test\target1-noMods.csv";
            outputFolderBox.Text = @"C:\Users\Anna\Desktop\NeuQuant CSMSL Test";
            rtWindow.Value = 1;
            channels.Value = 2;
            signalToNoiseThreshold.Value = 3;
            HILAC.Checked = true;
            noiseBandCap.Checked = false;
            coalescence.Checked = false;
        }

        private void Run()
        {
            RTWINDOW = (double) rtWindow.Value;
            NUMCHANNELS = (int) channels.Value;
            FIRSTSEARCHDONE = false;
            CHECKPAIRSPACING = true;
            NOISEBANDCAP = noiseBandCap.Checked;
            PEAKCOALESCENCE = coalescence.Checked;
            QUANTFILTER = true;
            CORRECTPROLINE = false;
            MULTIINJECT = false;
            NUMISOTOPES = 3;
            MINIMUMSN = (double) signalToNoiseThreshold.Value;
            CALCIUM = false;
            PEAKSNPAIRS = new Dictionary<int, List<double>>();

            if (!HILAC.Checked)
            {
                NEUCODE = false;
                if (NUMCHANNELS == 2)
                {
                    NUMISOTOPOLOGUES = 1;
                    NUMCLUSTERS = 2;
                    SILAC_DUPLEX_LYS = true;
                }
                else if (NUMCHANNELS == 3)
                {
                    NUMISOTOPOLOGUES = 2;
                    NUMCLUSTERS = 1;
                    NEUCODE_ARG_PROLINECONVERSION = true;
                    NUMCHANNELS = 2;
                }
            }
            else
            {
                NEUCODE = true;
                if (NUMCHANNELS == 2)
                {
                    NUMISOTOPOLOGUES = 2;
                    NUMCLUSTERS = 1;
                    NEUCODE_DUPLEX_LYS = true;
                    //NEUCODE_DUPLEX_LYS1 = true;
                }
                else if (NUMCHANNELS == 4)
                {
                    NUMISOTOPOLOGUES = 4;
                    NUMCLUSTERS = 1;
                    NEUCODE_4PLEX_LIGHT = true;
                    //NEUCODE_4PLEX_MEDIUM = true;
                    //NEUCODE_4PLEX_HEAVY = true;
                }
                else if (NUMCHANNELS == 6)
                {
                    NUMISOTOPOLOGUES = 2;
                    NUMCLUSTERS = 3;
                    NEUCODE_SIXPLEX_MTRAQ = true;
                    //NEUCODE_SIXPLEX_ARG = true;
                }
                else if (NUMCHANNELS == 12)
                {
                    NUMISOTOPOLOGUES = 4;
                    NUMCLUSTERS = 6;
                    NEUCODE_12PLEX = true;
                }
            }

            Console.WriteLine("starting");

            //Cycle through .csv files to make a list of identified peptides and properties
            Dictionary<string, PeptideID> allPeptides = new Dictionary<string, PeptideID>();
            List<PrecursorPPM> PRECURSORPPM = new List<PrecursorPPM>();
            RAWFILES = new Dictionary<string, MSDataFile>();
            List<CoalescenceCheck> INTENSITY_MISSINGCHANNEL = new List<CoalescenceCheck>();
            List<Spacing> spacings = new List<Spacing>();
            readCsvInputFile(allPeptides);         
            MSDataScan currentFullScan;
            int rawFileCount = 0;
            int totalRawFiles = RAWFILES.Count;

            Console.WriteLine("calculating systematic error");

            //Calculate systematic error by looking for monoisotopes in MS scan preceding best MS/MS scan
            foreach (MSDataFile rawFile in RAWFILES.Values)
            {
                rawFileCount++;
                rawFile.Open();
                Console.WriteLine("calculating ppm error in raw file " + rawFileCount + " of " + totalRawFiles);
                
                foreach (PeptideID peptide in allPeptides.Values)
                {
                    if (peptide.numLabels > 0)
                    {
                        List<PeptideSpectralMatch> psms = null;
                        if (peptide.PSMs.TryGetValue(rawFile, out psms))
                        {
                            // Sort PSMs by E-value
                            psms.OrderBy(psm => psm.EValue);
                            peptide.precursorPPMError(rawFile, PRECURSORPPM);
                        }
                    }
                }
            }
            
            PRECURSORPPM.Sort();
            
            // PPM error outputs can be printed out to a file to assess error distributions
            /*StreamWriter ppmWriter = new StreamWriter(Path.Combine(outputFolderBox.Text, Path.GetFileNameWithoutExtension(csvInputBox.Text) + "_ppm.csv"));
            string header = ("Peptide, Charge, PPM Error, E-value");
            ppmWriter.WriteLine(header);
            Console.WriteLine("writing output");

            foreach (PrecursorPPM ppm in PRECURSORPPM)
            {
                ppmWriter.WriteLine("{0},{1},{2},{3}", ppm.Peptide, ppm.Charge, ppm.Ppm, ppm.EValue);
            }

            ppmWriter.Close();*/

            SYSTEMATICERROR = PRECURSORPPM.ElementAt(PRECURSORPPM.Count / 2).Ppm; //Set systematic error as the median value of all light and heavy precursors
            
            Console.WriteLine("systematic error: " + SYSTEMATICERROR);
            PRECURSORPPM.Clear();

            Console.WriteLine("applying systematic error");
            FIRSTSEARCHDONE = true;

            Console.WriteLine("searching RAW file for pairs");
            rawFileCount = 0;
            foreach (MSDataFile rawFile in RAWFILES.Values)
            {
                rawFileCount++;
                rawFile.Open();
                Console.WriteLine("finding pairs in raw file " + rawFileCount + " of " + totalRawFiles);

                foreach (PeptideID uniquePeptide in allPeptides.Values)
                {
                    // Only consider peptides that contain at least one label
                    if (uniquePeptide.numLabels > 0)
                    {
                        // Only consider peptides that were identified in this raw file
                        List<PeptideSpectralMatch> psms = null;
                        if (uniquePeptide.PSMs.TryGetValue(rawFile, out psms))
                        {
                            // Find the appropriate precursor scan(s)
                            if (MULTIINJECT)
                            {
                                if (psms != null && psms.Count > 0)
                                {
                                    foreach (PeptideSpectralMatch psm in psms)
                                    {
                                        int correctScan = -1;
                                        bool stop = false;
                                        int scanCount = psm.ScanNumber + 1;
                                        while (scanCount < rawFile.LastSpectrumNumber && !stop)
                                        {
                                            if (rawFile.GetMsScan(scanCount).MsnOrder == 2 && ((MsnDataScan)rawFile.GetMsScan(scanCount)).DissociationType == DissociationType.HCD)
                                            {
                                                correctScan = scanCount;
                                                stop = true;
                                            }
                                            scanCount++;
                                        }
                                        currentFullScan = rawFile.GetMsScan(correctScan);
                                        uniquePeptide.findPeaks(currentFullScan, rawFile, psm.Charge);
                                    }
                                    if (CHECKPAIRSPACING)
                                    {
                                        uniquePeptide.checkPairSpacing(rawFile, spacings);
                                    }
                                }
                            }
                            else
                            {
                                int charge = psms[0].Charge;
                                uniquePeptide.calculateScanRange(rawFile, RTWINDOW);
                                if (uniquePeptide.fullScanList != null)
                                {
                                    uniquePeptide.fullScanList.Sort();
                                    foreach (int scanNumber in uniquePeptide.fullScanList)
                                    {
                                        currentFullScan = rawFile[scanNumber];
                                        uniquePeptide.findPeaks(currentFullScan, rawFile, charge);
                                    }

                                    if (CHECKPAIRSPACING)
                                    {
                                        uniquePeptide.checkPairSpacing(rawFile, spacings);
                                    }
                                }
                            }
                        }
                    }
                } 
                //rawFile.Close();
            }

            // Pair spacing outputs can be printed out to a file to assess spacing distributions
            /*StreamWriter spacingWriter = new StreamWriter(Path.Combine(outputFolderBox.Text, Path.GetFileNameWithoutExtension(csvInputBox.Text) + "_spacings.csv"));
            string header1 = ("Theo Spacing, Exp Spacing, Charge, MZ");
            spacingWriter.WriteLine(header1);
            Console.WriteLine("writing output");

            foreach (Spacing spacing in spacings)
            {
                spacingWriter.WriteLine("{0},{1},{2},{3}", spacing.theoSpacing, spacing.spacing, spacing.charge, spacing.MZ);
            }

            spacingWriter.Close();*/

            // Coalescence outputs can be printed out to a file to assess the impact of precursor signal-to-noise on observed coalescence
            /*StreamWriter coal = new StreamWriter(Path.Combine(outputFolderBox.Text, Path.GetFileNameWithoutExtension(csvInputBox.Text) + "_coalescence.csv"));
            string header1 = ("# Peaks, S/N");
            coal.WriteLine(header1);
            foreach (int numPeaks in PEAKSNPAIRS.Keys)
            {
                foreach (double signalNoise in PEAKSNPAIRS[numPeaks])
                {
                    coal.WriteLine("{0},{1}", numPeaks, signalNoise);
                }
            }

            coal.Close();*/
            
            if (NOISEBANDCAP && PEAKCOALESCENCE)
            {
                Console.WriteLine("calculating coalescence threshold");
                MAXIMUMDNL = calculateCoalescenceThreshold(allPeptides, 0.1, INTENSITY_MISSINGCHANNEL);
                Console.WriteLine("intensity threshold: " + MAXIMUMDNL);
            }

            // Validate NeuCode pairs in which one channel was found to be missing -- apply noise level or discard due to coalescence
            rawFileCount = 0;
            foreach (MSDataFile rawFile in RAWFILES.Values)
            {
                rawFileCount++;
                rawFile.Open();
                Console.WriteLine("quantifying pairs in raw file " + rawFileCount + " of " + totalRawFiles);

                foreach (PeptideID uniquePeptide in allPeptides.Values)
                {
                    if (uniquePeptide.numLabels > 0)
                    {
                        List<PeptideSpectralMatch> psms = null;
                        if (uniquePeptide.PSMs.TryGetValue(rawFile, out psms))
                        {
                            if (PEAKCOALESCENCE)
                            {                                
                                uniquePeptide.checkPairCoalescence(rawFile);
                            }

                            if (NOISEBANDCAP)
                            {                                
                                uniquePeptide.applyNoise(rawFile);
                            }
                        }
                    }
                }
            }

            foreach (PeptideID uniquePeptide in allPeptides.Values)
            {
                if (uniquePeptide.numLabels > 0)
                {
                    uniquePeptide.quantify();
                }
                else
                {
                    uniquePeptide.finalQuantified = new int[NUMCLUSTERS];
                    uniquePeptide.quantifiedNoiseIncluded = new bool[NUMCLUSTERS];
                    
                    for (int c = 0; c < NUMCLUSTERS; c++)
                    {
                        uniquePeptide.finalQuantified[c] = 0;
                        uniquePeptide.quantifiedNoiseIncluded[c] = false;
                    }
                }
            }

            Console.WriteLine("writing output file");
            writeCsvOutputFile(allPeptides);
        }

        private void writeCsvOutputFile(Dictionary<string, PeptideID> allPeptides)
        {
            string outputName;
            if (NEUCODE)
            {
                outputName = Path.Combine(outputFolderBox.Text, Path.GetFileNameWithoutExtension(csvInputBox.Text) + "_NeuCode_Quant.csv");
            }
            else
            {
                outputName = Path.Combine(outputFolderBox.Text, Path.GetFileNameWithoutExtension(csvInputBox.Text) + "_SILAC_Quant.csv");
            }
            StreamWriter writer1 = new StreamWriter(outputName);

            if (NUMCHANNELS == 2)
            {
                string header1 = ("Scan number(s), Raw File(s), Charge State(s), Best Scan Number, Peptide, # Labels, Theo Mass 1, Theo Mass 2, Adjusted Mass 1, Adjusted Mass 2, Quantified MS1 Scans, # Total Measurements, Total Intensity 1, Total Intensity 2, Quantified no NBC MS1 Scans, # no NBC Measurements, No NBC Intensity 1, No NBC Intensity 2, Coalescence Detected?, Ratio H/L, Ratio H/L Count, Missing Channels Quantified");
                writer1.WriteLine(header1);

                foreach (PeptideID peptide in allPeptides.Values)
                {
                    writer1.WriteLine("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21}",
                        peptide.scanNumbers, peptide.rawFiles, peptide.chargeStates, peptide.bestPSMScanNumber, peptide.sequence, peptide.numLabels,
                        peptide.theoMasses[0, 0], peptide.theoMasses[1, 0],
                        peptide.adjustedTheoMasses[0, 0], peptide.adjustedTheoMasses[1, 0],
                        peptide.countAllPairs, peptide.countAllIsotopes[0], peptide.totalIntensity[0, NUMISOTOPES], peptide.totalIntensity[1, NUMISOTOPES],
                        peptide.countCompletePairs, peptide.countCompleteIsotopes[0], peptide.completeTotalIntensity[0, NUMISOTOPES], peptide.completeTotalIntensity[1, NUMISOTOPES],
                        peptide.coalescenceDetected, peptide.heavyToLightRatioSum[0, 0], peptide.finalQuantified[0], peptide.quantifiedNoiseIncluded[0]);
                }
            }
            else if (NUMCHANNELS == 4)
            {
                string header1 = ("Scan number(s), Raw File(s), Charge State(s), Best Scan Number, Peptide, # Labels, Theo Mass 1, Theo Mass 2, Theo Mass 3, Theo Mass 4, Adjusted Mass 1, Adjusted Mass 2, Adjusted Mass 3, Adjusted Mass 4, Quantified MS1 Scans, # Total Measurements, Total Intensity 1, Total Intensity 2, Total Intensity 3, Total Intensity 4, Quantified no NBC MS1 Scans, # no NBC Measurements, No NBC Intensity 1, No NBC Intensity 2, No NBC Intensity 3, No NBC Intensity 4, Coalescence Detected?, Ratio 1, Ratio 2, Ratio 3, Ratio Count, Missing Channels Quantified");
                writer1.WriteLine(header1);

                foreach (PeptideID peptide in allPeptides.Values)
                {
                    writer1.WriteLine("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23},{24},{25},{26},{27},{28},{29},{30},{31}",
                        peptide.scanNumbers, peptide.rawFiles, peptide.chargeStates, peptide.bestPSMScanNumber, peptide.sequence, peptide.numLabels,
                        peptide.theoMasses[0, 0], peptide.theoMasses[1, 0], peptide.theoMasses[2,0], peptide.theoMasses[3,0],
                        peptide.adjustedTheoMasses[0, 0], peptide.adjustedTheoMasses[1, 0], peptide.adjustedTheoMasses[2,0], peptide.adjustedTheoMasses[3,0],
                        peptide.countAllPairs, peptide.countAllIsotopes[0], peptide.totalIntensity[0, NUMISOTOPES], peptide.totalIntensity[1, NUMISOTOPES], peptide.totalIntensity[2,NUMISOTOPES], peptide.totalIntensity[3,NUMISOTOPES],
                        peptide.countCompletePairs, peptide.countCompleteIsotopes[0], peptide.completeTotalIntensity[0, NUMISOTOPES], peptide.completeTotalIntensity[1, NUMISOTOPES], peptide.completeTotalIntensity[2,NUMISOTOPES], peptide.completeTotalIntensity[3,NUMISOTOPES],
                        peptide.coalescenceDetected, peptide.heavyToLightRatioSum[0,0], peptide.heavyToLightRatioSum[1,0], peptide.heavyToLightRatioSum[2,0], peptide.finalQuantified[0], peptide.quantifiedNoiseIncluded[0]);
                }
            }
            else if (NEUCODE_12PLEX)
            {
                string header2 = ("Scan number(s), Raw File(s), Charge State(s), Peptide, # Labels, Theo Mass 1, Theo Mass 2, Theo Mass 3, Theo Mass 4, Theo Mass 5, Theo Mass 6, Theo Mass 7, Theo Mass 8, Theo Mass 9, Theo Mass 10, Theo Mass 11, Theo Mass 12, Adjusted Mass 1, Adjusted Mass 2, Adjusted Mass 3, Adjusted Mass 4, Adjusted Mass 5, Adjusted Mass 6, Adjusted Mass 7, Adjusted Mass 8, Adjusted Mass 9, Adjusted Mass 10, Adjusted Mass 11, Adjusted Mass 12, Quantified MS1 Scans, # Total Measurements (Rep 1), # Total Measurements (Rep 2), # Total Measurements (Rep 3), Total Intensity 1, Total Intensity 2, Total Intensity 3, Total Intensity 4, Total Intensity 5, Total Intensity 6, Total Intensity 7, Total Intensity 8, Total Intensity 9, Total Intensity 10, Total Intensity 11, Total Intensity 12, Quantified no NBC MS1 Scans, # No NBC Measurements (Rep 1), # No NBC Measurements (Rep 2), # No NBC Measurements (Rep 3), No NBC Intensity 1, No NBC Intensity 2, No NBC Intensity 3, No NBC Intensity 4, No NBC Intensity 5, No NBC Intensity 6, No NBC Intensity 7, No NBC Intensity 8, No NBC Intensity 9, No NBC Intensity 10, No NBC Intensity 11, No NBC Intensity 12, Coalescence Detected?, Ratio 1 (Rep 1), Ratio 2 (Rep 1), Ratio 3 (Rep 1), Ratio Count (Rep 1), Missing Channels Quantified (Rep 1), Ratio 1 (Rep 2), Ratio 2 (Rep 2), Ratio 3 (Rep 2), Ratio Count (Rep 2), Missing Channels Quantified (Rep 2), Ratio 1 (Rep 3), Ratio 2 (Rep 3), Ratio 3 (Rep 3), Ratio Count (Rep 3), Missing Channels Quantified (Rep 3)");
                writer1.WriteLine(header2);

                foreach (PeptideID peptide in allPeptides.Values)
                {
                    writer1.WriteLine("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23},{24},{25},{26},{27},{28},{29},{30},{31},{32},{33},{34},{35},{36},{37},{38},{39},{40},{41},{42},{43},{44},{45},{46},{47},{48},{49},{50},{51},{52},{53},{54},{55},{56},{57},{58},{59},{60},{61},{62},{63},{64},{65},{66},{67},{68},{69},{70},{71},{72},{73},{74},{75},{76}",
                        peptide.scanNumbers, peptide.rawFiles, peptide.chargeStates, peptide.sequence, peptide.numLabels,
                        peptide.theoMasses[0, 0], peptide.theoMasses[1, 0], peptide.theoMasses[2, 0], peptide.theoMasses[3, 0], peptide.theoMasses[4, 0], peptide.theoMasses[5, 0], peptide.theoMasses[6, 0], peptide.theoMasses[7, 0], peptide.theoMasses[8, 0], peptide.theoMasses[9, 0], peptide.theoMasses[10, 0], peptide.theoMasses[11, 0],
                        peptide.adjustedTheoMasses[0, 0], peptide.adjustedTheoMasses[1, 0], peptide.adjustedTheoMasses[2, 0], peptide.adjustedTheoMasses[3, 0], peptide.adjustedTheoMasses[4, 0], peptide.adjustedTheoMasses[5, 0], peptide.adjustedTheoMasses[6, 0], peptide.adjustedTheoMasses[7, 0], peptide.adjustedTheoMasses[8, 0], peptide.adjustedTheoMasses[9, 0], peptide.adjustedTheoMasses[10, 0], peptide.adjustedTheoMasses[11, 0],
                        peptide.countAllPairs, peptide.countAllIsotopes[0], peptide.countAllIsotopes[1], peptide.countAllIsotopes[2],
                        peptide.totalIntensity[0, NUMISOTOPES], peptide.totalIntensity[1, NUMISOTOPES], peptide.totalIntensity[2, NUMISOTOPES], peptide.totalIntensity[3, NUMISOTOPES], peptide.totalIntensity[4, NUMISOTOPES], peptide.totalIntensity[5, NUMISOTOPES], peptide.totalIntensity[6, NUMISOTOPES], peptide.totalIntensity[7, NUMISOTOPES], peptide.totalIntensity[8, NUMISOTOPES], peptide.totalIntensity[9, NUMISOTOPES], peptide.totalIntensity[10, NUMISOTOPES], peptide.totalIntensity[11, NUMISOTOPES],
                        peptide.countCompletePairs, peptide.countCompleteIsotopes[0], peptide.countCompleteIsotopes[1], peptide.countCompleteIsotopes[2],
                        peptide.completeTotalIntensity[0, NUMISOTOPES], peptide.completeTotalIntensity[1, NUMISOTOPES], peptide.completeTotalIntensity[2, NUMISOTOPES], peptide.completeTotalIntensity[3, NUMISOTOPES], peptide.completeTotalIntensity[4, NUMISOTOPES], peptide.completeTotalIntensity[5, NUMISOTOPES], peptide.completeTotalIntensity[6, NUMISOTOPES], peptide.completeTotalIntensity[7, NUMISOTOPES], peptide.completeTotalIntensity[8, NUMISOTOPES], peptide.completeTotalIntensity[9, NUMISOTOPES], peptide.completeTotalIntensity[10, NUMISOTOPES], peptide.completeTotalIntensity[11, NUMISOTOPES],
                        peptide.coalescenceDetected, 
                        peptide.heavyToLightRatioSum[0,0], peptide.heavyToLightRatioSum[1,0], peptide.heavyToLightRatioSum[2,0], peptide.finalQuantified[0], peptide.quantifiedNoiseIncluded[0],
                        peptide.heavyToLightRatioSum[3,0], peptide.heavyToLightRatioSum[4,0], peptide.heavyToLightRatioSum[5,0], peptide.finalQuantified[1], peptide.quantifiedNoiseIncluded[1],
                        peptide.heavyToLightRatioSum[6,0], peptide.heavyToLightRatioSum[7,0], peptide.heavyToLightRatioSum[8,0], peptide.finalQuantified[2], peptide.quantifiedNoiseIncluded[2]);
                }
            }
            /*else if (NEUCODE_ARG_PROLINECONVERSION)
            {
                string header2 = ("Scan number(s), Raw File(s), Charge State(s), Best Scan Number, Peptide, # K, # R, No Pro Light Theo Mass, No Pro Heavy Theo Mass, Pro Light Theo Mass, Pro Heavy Theo Mass, No Pro Adjusted Light Theo Mass, No Pro Adjusted Heavy Theo Mass, Pro Adjusted Light Theo Mass, Pro Adjusted Heavy Theo Mass, Quantified MS1 Scans, # Total Measurements, No Pro Light Intensity, No Pro Heavy Intensity, Pro Light Intensity, Pro Heavy Intensity, Quantified no NBC MS1 Scans, # no NBC Measurements, No Pro Light no NBC Intensity, No Pro Heavy no NBC Intensity, Pro Light no NBC Intensity, Pro Heavy no NBC Intensity, Coalescence Detected?, No Pro Ratio H/L, Pro Ratio H/L, Ratio H/L Count, Missing Channels Quantified?");
                writer1.WriteLine(header2);

                foreach (PeptideID peptide in allPeptides.Values)
                {
                    writer1.WriteLine("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23},{24},{25},{26},{27},{28},{29},{30},{31}",
                        peptide.scanNumbers, peptide.rawFiles, peptide.chargeStates, peptide.bestPSMScanNumber, peptide.sequence, peptide.lysineCount, peptide.arginineCount,
                        peptide.theoMasses[0, 0], peptide.theoMasses[1, 0], peptide.theoMasses[2, 0], peptide.theoMasses[3, 0],
                        peptide.adjustedTheoMasses[0, 0], peptide.adjustedTheoMasses[1, 0], peptide.adjustedTheoMasses[2, 0], peptide.adjustedTheoMasses[3, 0],
                        peptide.countAllPairs, peptide.countAllIsotopes, peptide.totalIntensity[0, NUMISOTOPES], peptide.totalIntensity[1, NUMISOTOPES], peptide.totalIntensity[2, NUMISOTOPES], peptide.totalIntensity[3, NUMISOTOPES],
                        peptide.countNoNBCPairs, peptide.countNoNBCIsotopes, peptide.noNBCTotalIntensity[0, NUMISOTOPES], peptide.noNBCTotalIntensity[1, NUMISOTOPES], peptide.noNBCTotalIntensity[2, NUMISOTOPES], peptide.noNBCTotalIntensity[3, NUMISOTOPES],
                        peptide.coalescenceDetected, peptide.heavyToLightRatioSum[0, 0], peptide.heavyToLightRatioSum[1, 0], peptide.finalQuantified, peptide.quantifiedNBC);
                }
            }*/
            else if (NEUCODE_SIXPLEX_MTRAQ)
            {
                string header2 = ("Scan number(s), Raw File(s), Charge State(s), Best Scan Number, Peptide, # Labels, Theo Mass 1, Theo Mass 2, Theo Mass 3, Theo Mass 4, Theo Mass 5, Theo Mass 6, Adjusted Mass 1, Adjusted Mass 2, Adjusted Mass 3, Adjusted Mass 4, Adjusted Mass 5, Adjusted Mass 6, Quantified MS1 Scans, # Total Measurements (Rep 1), # Total Measurements (Rep 2), # Total Measurements (Rep 3), Total Intensity 1, Total Intensity 2, Total Intensity 3, Total Intensity 4, Total Intensity 5, Total Intensity 6, Quantified no NBC MS1 Scans, # No NBC Measurements (Rep 1), # No NBC Measurements (Rep 2), # No NBC Measurements (Rep 3), No NBC Intensity 1, No NBC Intensity 2, No NBC Intensity 3, No NBC Intensity 4, No NBC Intensity 5, No NBC Intensity 6, Coalescence Detected?, Ratio (Rep 1), Ratio Count (Rep 1), Missing Channels Quantified (Rep 1), Ratio (Rep 2), Ratio Count (Rep 2), Missing Channels Quantified (Rep 2), Ratio (Rep 3), Ratio Count (Rep 3), Missing Channels Quantified (Rep 3)");
                writer1.WriteLine(header2);

                foreach (PeptideID peptide in allPeptides.Values)
                {
                    writer1.WriteLine("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23},{24},{25},{26},{27},{28},{29},{30},{31},{32},{33},{34},{35},{36},{37},{38},{39},{40},{41},{42},{43},{44},{45},{46},{47}",
                        peptide.scanNumbers, peptide.rawFiles, peptide.chargeStates, peptide.bestPSMScanNumber, peptide.sequence, peptide.numLabels,
                        peptide.theoMasses[0, 0], peptide.theoMasses[1, 0], peptide.theoMasses[2, 0], peptide.theoMasses[3, 0], peptide.theoMasses[4, 0], peptide.theoMasses[5, 0],
                        peptide.adjustedTheoMasses[0, 0], peptide.adjustedTheoMasses[1, 0], peptide.adjustedTheoMasses[2, 0], peptide.adjustedTheoMasses[3, 0], peptide.adjustedTheoMasses[4, 0], peptide.adjustedTheoMasses[5, 0],
                        peptide.countAllPairs, peptide.countAllIsotopes[0], peptide.countAllIsotopes[1], peptide.countAllIsotopes[2],
                        peptide.totalIntensity[0, NUMISOTOPES], peptide.totalIntensity[1, NUMISOTOPES], peptide.totalIntensity[2, NUMISOTOPES], peptide.totalIntensity[3, NUMISOTOPES], peptide.totalIntensity[4, NUMISOTOPES], peptide.totalIntensity[5, NUMISOTOPES],
                        peptide.countCompletePairs, peptide.countCompleteIsotopes[0], peptide.countCompleteIsotopes[1], peptide.countCompleteIsotopes[2],
                        peptide.completeTotalIntensity[0, NUMISOTOPES], peptide.completeTotalIntensity[1, NUMISOTOPES], peptide.completeTotalIntensity[2, NUMISOTOPES], peptide.completeTotalIntensity[3, NUMISOTOPES], peptide.completeTotalIntensity[4, NUMISOTOPES], peptide.completeTotalIntensity[5, NUMISOTOPES],
                        peptide.coalescenceDetected, 
                        peptide.heavyToLightRatioSum[0, 0], peptide.finalQuantified[0], peptide.quantifiedNoiseIncluded[0],
                        peptide.heavyToLightRatioSum[1, 0], peptide.finalQuantified[1], peptide.quantifiedNoiseIncluded[1],
                        peptide.heavyToLightRatioSum[2, 0], peptide.finalQuantified[2], peptide.quantifiedNoiseIncluded[2]);
                }
            }
            else
            {
                string header3 = ("Scan number(s), Raw File(s), Charge State(s), Best Scan Number, Peptide, # K, # R, Light R Light Theo Mass, Light R Heavy Theo Mass, Light R Adjusted Light Theo Mass, Light R Adjusted Heavy Theo Mass, Medium R Light Theo Mass, Medium R Heavy Theo Mass, Medium R Adjusted Light Theo Mass, Medium R Adjusted Heavy Theo Mass, Heavy R Light Theo Mass, Heavy R Heavy Theo Mass, Heavy R Adjusted Light Theo Mass, Heavy R Adjusted Heavy Theo Mass, Quantified MS1 Scans, # Total Measurements, Total Light R Light Intensity, Total Light R Heavy Intensity, Total Medium R Light Intensity, Total Medium R Heavy Intensity, Total Heavy R Light Intensity, Total Heavy R Heavy Intensity, Quantified no NBC MS1 Scans, # No NBC Measurements (Rep 1), # No NBC Measurements (Rep 2), # No NBC Measurements (Rep 3), Light R Light no NBC Intensity, Light R Heavy no NBC Intensity, Medium R Light no NBC Intensity, Medium R Heavy no NBC Intensity, Heavy R Light no NBC Intensity, Coalescence Detected?, Ratio (Rep 1), Ratio Count (Rep 1), Missing Channels Quantified (Rep 1), Ratio (Rep 2), Ratio Count (Rep 2), Missing Channels Quantified (Rep 2), Ratio (Rep 3), Ratio Count (Rep 3), Missing Channels Quantified (Rep 3)");
                writer1.WriteLine(header3);

                foreach (PeptideID peptide in allPeptides.Values)
                {
                    if (peptide.arginineCount == 0)
                    {
                        writer1.WriteLine("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23},{24},{25},{26},{27},{28},{29},{30},{31},{32},{33},{34},{35},{36},{37},{38},{39},{40}",
                        peptide.scanNumbers, peptide.rawFiles, peptide.chargeStates, peptide.bestPSMScanNumber, peptide.sequence, peptide.lysineCount, peptide.arginineCount,
                        peptide.theoMasses[0, 0], peptide.theoMasses[1, 0], 0, 0, 0, 0,
                        peptide.adjustedTheoMasses[0, 0], peptide.adjustedTheoMasses[1, 0], 0, 0, 0, 0,
                        peptide.countAllPairs, peptide.countAllIsotopes, peptide.totalIntensity[0, NUMISOTOPES], peptide.totalIntensity[1, NUMISOTOPES], 0, 0, 0, 0,
                        peptide.countCompletePairs, peptide.countCompleteIsotopes, peptide.completeTotalIntensity[0, NUMISOTOPES], peptide.completeTotalIntensity[1, NUMISOTOPES], 0, 0, 0, 0,
                        peptide.coalescenceDetected, peptide.heavyToLightRatioSum[0, 0], 0, 0, peptide.finalQuantified, peptide.quantifiedNoiseIncluded);
                    }
                    else
                    {
                        writer1.WriteLine("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23},{24},{25},{26},{27},{28},{29},{30},{31},{32},{33},{34},{35},{36},{37},{38},{39},{40}",
                        peptide.scanNumbers, peptide.rawFiles, peptide.chargeStates, peptide.bestPSMScanNumber, peptide.sequence, peptide.lysineCount, peptide.arginineCount,
                        peptide.theoMasses[0, 0], peptide.theoMasses[1, 0], peptide.theoMasses[2, 0], peptide.theoMasses[3, 0], peptide.theoMasses[4, 0], peptide.theoMasses[5, 0],
                        peptide.adjustedTheoMasses[0, 0], peptide.adjustedTheoMasses[1, 0], peptide.adjustedTheoMasses[2, 0], peptide.adjustedTheoMasses[3, 0], peptide.adjustedTheoMasses[4, 0], peptide.adjustedTheoMasses[5, 0],
                        peptide.countAllPairs, peptide.countAllIsotopes, peptide.totalIntensity[0, NUMISOTOPES], peptide.totalIntensity[1, NUMISOTOPES], peptide.totalIntensity[2, NUMISOTOPES], peptide.totalIntensity[3, NUMISOTOPES], peptide.totalIntensity[4, NUMISOTOPES], peptide.totalIntensity[5, NUMISOTOPES],
                        peptide.countCompletePairs, peptide.countCompleteIsotopes, peptide.completeTotalIntensity[0, NUMISOTOPES], peptide.completeTotalIntensity[1, NUMISOTOPES], peptide.completeTotalIntensity[2, NUMISOTOPES], peptide.completeTotalIntensity[3, NUMISOTOPES], peptide.completeTotalIntensity[4, NUMISOTOPES], peptide.completeTotalIntensity[5, NUMISOTOPES],
                        peptide.coalescenceDetected, peptide.heavyToLightRatioSum[0, 0], peptide.heavyToLightRatioSum[1, 0], peptide.heavyToLightRatioSum[2, 0], peptide.finalQuantified, peptide.quantifiedNoiseIncluded);
                    }
                }
            }
            writer1.Close();
        }

        /* Reads in peptide information from a CSV file
         * Necessary information: spectrum number, charge, peptide sequence, E-value, raw file
         */
        private void readCsvInputFile(Dictionary<string, PeptideID> allPeptides)
        {
            if (NEUCODE)
            {
                if (NUMISOTOPOLOGUES != 2 && NUMISOTOPOLOGUES != 4)                             
                {
                    Console.WriteLine("invalid number of multiplexing channels for selected NeuCode approach");
                    Application.Exit();
                }
            }
            else
            {
                if (NUMCHANNELS != 2)
                {
                    Console.WriteLine("invalid number of multiplexing channels for selected SILAC approach");
                    Application.Exit();
                }
            }

            //Amino acids

            NamedChemicalFormula K13C15N = NamedChemicalFormula.AddModification("C-6C{13}6N-2N{15}2", "Lys +8 13C6 15N2");
            NamedChemicalFormula K2H = NamedChemicalFormula.AddModification("H-8H{2}8", "Lys +8 2H8");
            NamedChemicalFormula R13C = NamedChemicalFormula.AddModification("C-6C{13}6", "Arg +6 13C6");
            NamedChemicalFormula R13C15N = NamedChemicalFormula.AddModification("C-6C{13}6N-4N{15}2", "Arg + 10 13C615N4");
            NamedChemicalFormula K13C = NamedChemicalFormula.AddModification("C-1C{13}", "Lys +1 13C");
            NamedChemicalFormula K15N = NamedChemicalFormula.AddModification("N-1N{15}", "Lys +1 15N");

            //Chemical labels
            NamedChemicalFormula lightmTRAQ = NamedChemicalFormula.AddModification("H{1}12C{12}7N{14}2O{16}1", "mTRAQ L");
            NamedChemicalFormula lightmTRAQK13C15N = NamedChemicalFormula.AddModification("H{1}12C{12}7N{14}2O{16}1C-6C{13}6N-2N{15}2", "mTRAQ L Lys +8 13C6 15N2");
            NamedChemicalFormula lightmTRAQK2H = NamedChemicalFormula.AddModification("H{1}12C{12}7N{14}2O{16}1H-8H{2}8", "mTRAQ L Lys +8 2H8");
            NamedChemicalFormula mediummTRAQ = NamedChemicalFormula.AddModification("H{1}12C{12}4C{13}3N{14}1N{15}1O{16}1", "mTRAQ M");
            NamedChemicalFormula mediummTRAQK13C15N = NamedChemicalFormula.AddModification("H{1}12C{12}4C{13}3N{14}1N{15}1O{16}1C-6C{13}6N-2N{15}2", "mTRAQ M Lys +8 13C6 15N2");
            NamedChemicalFormula mediummTRAQK2H = NamedChemicalFormula.AddModification("H{1}12C{12}4C{13}3N{14}1N{15}1O{16}1H-8H{2}8", "mTRAQ M Lys +8 2H8");
            NamedChemicalFormula heavymTRAQ = NamedChemicalFormula.AddModification("H{1}12C{12}1C{13}6N{15}2O{16}1", "mTRAQ H");
            NamedChemicalFormula heavymTRAQK13C15N = NamedChemicalFormula.AddModification("H{1}12C{12}1C{13}6N{15}2O{16}1C-6C{13}6N-2N{15}2", "mTRAQ H Lys +8 13C6 15N2");
            NamedChemicalFormula heavymTRAQK2H = NamedChemicalFormula.AddModification("H{1}12C{12}1C{13}6N{15}2O{16}1H-8H{2}8", "mTRAQ H Lys +8 2H8");
            NamedChemicalFormula lightASH1 = NamedChemicalFormula.AddModification("C{12}18H{1}31N{14}1O{16}5N{15}6", "4plex L1");
            NamedChemicalFormula lightASH2 = NamedChemicalFormula.AddModification("C{12}16H{1}31N{14}3O{16}5N{15}4C{13}2", "4plex L2");
            NamedChemicalFormula lightASH3 = NamedChemicalFormula.AddModification("C{12}14H{1}31N{14}5O{16}5N{15}2C{13}4", "4plex L3");
            NamedChemicalFormula lightASH4 = NamedChemicalFormula.AddModification("C{12}12H{1}31N{14}7O{16}5C{13}6", "4plex L4");
            NamedChemicalFormula mediumASH1 = NamedChemicalFormula.AddModification("C{12}14H{1}31N{14}1O{16}5N{15}6C{13}4", "4plex M1");
            NamedChemicalFormula mediumASH2 = NamedChemicalFormula.AddModification("C{12}12H{1}31N{14}1O{16}5N{15}4C{13}6", "4plex M2");
            NamedChemicalFormula mediumASH3 = NamedChemicalFormula.AddModification("C{12}10H{1}31N{14}1O{16}5N{15}2C{13}8", "4plex M3");
            NamedChemicalFormula mediumASH4 = NamedChemicalFormula.AddModification("C{12}8H{1}31N{14}1O{16}5C{13}10", "4plex M4");
            NamedChemicalFormula heavyASH1 = NamedChemicalFormula.AddModification("C{12}10H{1}31N{14}1O{16}5N{15}6C{13}8", "4plex H1");
            NamedChemicalFormula heavyASH2 = NamedChemicalFormula.AddModification("C{12}8H{1}31N{14}1O{16}5N{15}4C{13}10", "4plex H2");
            NamedChemicalFormula heavyASH3 = NamedChemicalFormula.AddModification("C{12}6H{1}31N{14}1O{16}5N{15}2C{13}12", "4plex H3");
            NamedChemicalFormula heavyASH4 = NamedChemicalFormula.AddModification("C{12}4H{1}31N{14}1O{16}5C{13}14", "4plex H4");
            NamedChemicalFormula lightCarbamyl = NamedChemicalFormula.AddModification("C{12}1N{15}1H{1}2O{16}1", "Carbamyl L");
            NamedChemicalFormula heavyCarbamyl = NamedChemicalFormula.AddModification("C{13}1N{14}1H{1}2O{16}1", "Carbamyl H");

            HashSet<string> rawFiles = new HashSet<string>();

            //Cycle through .csv file to make a list of identified peptides and properties
            CsvReader reader = new CsvReader(new StreamReader(csvInputBox.Text), true);
            using (reader)
            {
                PeptideID newPeptide;
                MSDataFile rawFile = null;
                Console.WriteLine("uploading peptide IDs");
                while (reader.ReadNextRecord())
                {
                    string basePathName = reader["Filename/id"].Substring(0, reader["Filename/id"].IndexOf("."));

                    if (RAWFILES.TryGetValue(basePathName, out rawFile))
                    {

                    }
                    else
                    {
                        rawFile = new ThermoRawFile(Path.Combine(rawFileBox.Text, basePathName + ".raw"));
                        RAWFILES.Add(basePathName, rawFile);
                    }                  

                    newPeptide = new PeptideID(int.Parse(reader["Spectrum number"]),
                                int.Parse(reader["Charge"]), double.Parse(reader["E-value"]),
                                reader["Peptide"], rawFile, reader["Mods"]);
                    checkAdd(allPeptides, newPeptide, (MSDataFile) rawFile, double.Parse(reader["E-value"]));                 
                }
                Console.WriteLine("done uploading peptide IDs: " + allPeptides.Count + " total unique sequences");
            }
        }

        /* Organizes the main peptide dictionary upon addition of each new peptide PSM
         */
        private void checkAdd(Dictionary<string, PeptideID> allPeptides, PeptideID peptide, MSDataFile rawFile, double eValue)
        {
            PeptideID peptideID = null;
            // Search main peptide dictionary
            if (allPeptides.TryGetValue(peptide.sequence, out peptideID))
            {
                // Peptide already in the dictionary
                List<PeptideSpectralMatch> PSMs = null;
                // Search peptide's PSM dictionary
                if (peptideID.PSMs.TryGetValue(rawFile, out PSMs))
                {
                    // Raw file already in dictionary -- update PSMs
                    PSMs.Add(peptide.PSM);
                }
                else 
                {
                    // Raw file not in PSM dictionary
                    peptideID.PSMs.Add(rawFile, new List<PeptideSpectralMatch>());
                    List<PeptideSpectralMatch> psms = null;
                    peptideID.PSMs.TryGetValue(rawFile, out psms);
                    psms.Add(peptide.PSM);
                    peptideID.PSMs[rawFile].Add(peptide.PSM);
                    peptideID.allPairs.Add(rawFile, new List<Pair>());
                    peptideID.completePairs.Add(rawFile, new List<Pair>());
                }

            }
            else
            {
                // Peptide not in main dictionary
                allPeptides.Add(peptide.sequence, peptide);
            }
        }

        /* Sort peptide max intensities into bins, calculate each bin's average missing channel frequency & set an intensity cutoff for coalescence
         */
        private double calculateCoalescenceThreshold(Dictionary<string, PeptideID> allPeptides, double binSize, List<CoalescenceCheck> list)
        {
            double intensityCutOff = 10000000.0;
            double intensityBinSize = binSize;
            CoalescenceCheck check;

            /*StreamWriter coalWriter = new StreamWriter(Path.Combine(outputFolderBox.Text, Path.GetFileNameWithoutExtension(csvInputBox.Text) + "_coalescencePlot.csv"));
            string header = ("Intensity Bin, # Peptides, Average Missing Channel Frequency");
            coalWriter.WriteLine(header);
            Console.WriteLine("writing output");*/
    
            //Retrieve each peptide's maximum intensity & missing channel frequency
            foreach (PeptideID uniquePeptide in allPeptides.Values)
            {
                double maximumIntensity = uniquePeptide.maximumIntensity;
                double missingChannelFrequency = uniquePeptide.missingChannelFrequency;
                
                if (uniquePeptide.maximumIntensity > 0 && missingChannelFrequency >= 0)
                {
                    check = new CoalescenceCheck(maximumIntensity, missingChannelFrequency);
                    list.Add(check);
                }
            }

            list.Sort();

            // Set up the intensity bins (log10)
            int intensityCount = list.Count;
            double minIntensity = list[0].intensity;
            double maxIntensity = list[intensityCount - 1].intensity;
            int numBins = (int)((maxIntensity - minIntensity) / intensityBinSize) + 1;

            List<double>[] intensityBins = new List<double>[numBins]; // Each array position represents a different intensity bin
            for (int i = 0; i < numBins; i++)
            {
                intensityBins[i] = new List<double>();
            }
            double[] binAverageFrequency = new double[numBins]; // Each array position represents the average missing channel frequency for the intensity bin

            double minIntensityRounded = System.Math.Round(minIntensity, 1);
            double firstBinMaximum = minIntensityRounded;
            double lastBinMaximum = firstBinMaximum + (double)(numBins * intensityBinSize);

            // Sort peptide missing channel frequencies into bins based on the peptide's maximum intensity
            double intensity;
            double frequency;
            int index;
            foreach (CoalescenceCheck coalCheck in list)
            {
                intensity = System.Math.Round(coalCheck.intensity, 2);
                frequency = coalCheck.missingChannelFrequency;
                if (intensity < firstBinMaximum)
                {
                    index = 0;
                }
                else if (intensity >= lastBinMaximum)
                {
                    index = numBins - 1;
                }
                else
                {
                    index = (int)((intensity - firstBinMaximum) / (intensityBinSize));
                }
                intensityBins[index].Add(frequency);
            }

            // Calculate the average missing channel frequency for each intensity bin & the minimum
            double averageFrequency;
            double sum;
            int count;
            List<double> freqs;
            double intensityValue;
            SortedList<double, double> sortedLowestFreqs = new SortedList<double, double>(); // Keep track of which 5 intensity bins produce the lowest frequencies
            sortedLowestFreqs.Add(1.0, 0.0);
            sortedLowestFreqs.Add(1.1, 0.0);
            sortedLowestFreqs.Add(1.2, 0.0);
            sortedLowestFreqs.Add(1.3, 0.0);
            sortedLowestFreqs.Add(1.4, 0.0);

            for (int i = 0; i < numBins; i++)
            {
                sum = 0;
                count = 0;
                freqs = intensityBins[i];

                // Calculate the bin average frequency
                foreach (double freq in freqs)
                {
                    sum += freq;
                    count++;
                }
                averageFrequency = sum / (double)count;
                intensityValue = (double)(i * intensityBinSize) + firstBinMaximum;

                //coalWriter.WriteLine("{0},{1},{2}", intensityValue, count, averageFrequency);
                
                // Check if bin average frequency is one of the 5 lowest out of all intensity bins
                if (averageFrequency < sortedLowestFreqs.ElementAt(4).Key && averageFrequency > 0)
                {
                    intensityValue = (double)(i * intensityBinSize) + firstBinMaximum;
                    sortedLowestFreqs.Add(averageFrequency, intensityValue);
                    sortedLowestFreqs.RemoveAt(5);
                }
                binAverageFrequency[i] = averageFrequency;
            }
            //coalWriter.Close();

            // Calculate and return the average intensity producing the minimum missing channel frequency
            sum = 0;
            foreach (KeyValuePair<double, double> pair in sortedLowestFreqs)
            {
                sum += pair.Value;
            }

            intensityCutOff = System.Math.Pow(10.0, sum / 5.0);

            //coalWriter.Close();
            
            return intensityCutOff;
        }

        private void RAWBrowse_Click(object sender, EventArgs e)
        {
            if (browseRawLocation.ShowDialog() == DialogResult.OK)
            {
                rawFileBox.Text = browseRawLocation.SelectedPath;
            }
        }

        private void CSVBrowse_Click(object sender, EventArgs e)
        {
            if (browseTargetInput.ShowDialog() == DialogResult.OK)
            {
                csvInputBox.Text = browseTargetInput.FileName;
            }
        }

        private void OutputBrowse_Click(object sender, EventArgs e)
        {
            if (browseOutputLocation.ShowDialog() == DialogResult.OK)
            {
                outputFolderBox.Text = browseOutputLocation.SelectedPath;
            }
        }

        private void start_Click(object sender, EventArgs e)
        {
            Run();
        }

        private void Form1_Load(object sender, EventArgs e)
        {

        }
    }
}
