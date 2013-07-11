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
        // Experiment constants
        public static double RTWINDOW;
        public static int NUMISOTOPES;
        public static double MINIMUMSN;
        public static MSDataFile RAWFILE;
        public static double SYSTEMATICERROR;
        public static bool FIRSTSEARCHDONE;
        public static Dictionary<string, MSDataFile> RAWFILES;
        public static Dictionary<string, Dictionary<int, Dictionary<Range<double>, double>>> INJECTIONTIMES;
        public static double MAXIMUMDNL;
        public static double THEORETICALSEPARATION;
        public static double QUANTRESOLUTION;
        public static double TOLERANCE;

        // Optional analyses
        public static bool CHECKPAIRSPACING;
        public static bool NOISEBANDCAP;
        public static bool PEAKCOALESCENCE;
        public static bool QUANTFILTER;
        public static bool CORRECTARGPROLINE;
        public static bool CORRECTLEUDLOSS;
        public static bool CORRECTLEUNLOSS;
        public static bool MULTIINJECT;
        public static bool AGCBINNING;
        //public static bool CALCIUM;
        public static bool CHECKPARTIALINCORPORATION;
        public static bool SEGMENTEDINJECTIONTIMES;
        //public static Dictionary<int, List<double>> PEAKSNPAIRS;
        //public static int QUANTCOUNT;
        
        // Experiment types
        public static bool NEUCODE;
        public static bool TRADITIONAL;
        public static bool ICAT;
        public static bool NEUCODE_DUPLEX_LYS8_36MDA;
        public static bool NEUCODE_DUPLEX_LEU7_18MDA;
        public static bool NEUCODE_TRIPLEX_LYS8_18MDA;
        public static bool NEUCODE_FOURPLEX_LYS8_12MDA;
        public static bool NEUCODE_SIXPLEX_LYS8_6MDA;
        public static bool NEUCODE_DUPLEX_LYS1_6MDA;
        public static bool NEUCODE_DUPLEX_CARBAMYL;
        public static bool SILAC_DUPLEX_LYSC;
        public static bool SILAC_DUPLEX_LYSCN;
        public static bool SILAC_DUPLEX_LYSH;
        public static bool SILAC_DUPLEX_LEUCN;
        public static bool SILAC_DUPLEX_LEUH;
        public static bool NEUCODE_SIXPLEX_MTRAQ;
        public static bool NEUCODE_SIXPLEX_ARG;
        public static bool NEUCODE_SIXPLEX_LEU;
        public static bool NEUCODE_12PLEX;
        public static bool NEUCODE_4PLEX_HEAVY;
        public static bool NEUCODE_4PLEX_MEDIUM;
        public static bool NEUCODE_4PLEX_LIGHT;
        public static bool NEUCODE_9PLEX_MTRAQ;
        public static bool NEUCODE_12PLEX_MTRAQ;
        public static bool NEUCODE_18PLEX_MTRAQ;
        public static int NUMCHANNELS;
        public static int NUMISOTOPOLOGUES;
        public static int NUMCLUSTERS;
        
        public Form1()
        {

            InitializeComponent();
            rawFileBox.Text = @"E:\Desktop\4-plex Yeast\Biological\RAW\4-plex";
            //listBox1.Items.Add(@"E:\Desktop\NeuCode 4-plex KGG\FDR\07092013 Josh B\DK06_95LDA_50ppm_trypsin_StaticK8_30k_search_peptides_OMSSAformat.csv");
            outputFolderBox.Text = @"E:\Desktop\4-plex Yeast\Biological\QUANT\4-plex";
            //csvInputBox.Text = @"E:\Desktop\NeuCode vs. TMT 4-plex\FDR\target-decoy\AEM_ITMS_CID_target.csv";
        }

        private void Run(string file)
        {
            RTWINDOW = (double) rtWindow.Value;
            NUMISOTOPES = (int)Isotopes.Value;
            MINIMUMSN = (double)signalToNoiseThreshold.Value;
            THEORETICALSEPARATION = (double)PeakSeparation.Value;
            QUANTRESOLUTION = (double)QuantResolution.Value;
            TOLERANCE = (double)searchTolerance.Value;
            CHECKPAIRSPACING = true;
            QUANTFILTER = true;
            NOISEBANDCAP = noiseBandCap.Checked;
            PEAKCOALESCENCE = coalescence.Checked;
            MULTIINJECT = MultipleInjections.Checked;
            AGCBINNING = AGCBins.Checked;

            Dictionary<string, PeptideID> allPeptides = new Dictionary<string, PeptideID>();
            List<PrecursorPPM> PRECURSORPPM = new List<PrecursorPPM>();
            RAWFILES = new Dictionary<string, MSDataFile>();
            List<CoalescenceCheck> INTENSITY_MISSINGCHANNEL = new List<CoalescenceCheck>();
            List<Spacing> spacings = new List<Spacing>();

            NEUCODE_DUPLEX_LYS8_36MDA = NeuCodeLys8Duplex.Checked && !mTRAQ.Checked && !Arg.Checked && !Leu.Checked;
            NEUCODE_DUPLEX_LEU7_18MDA = NeuCodeLeu7Duplex.Checked;
            NEUCODE_TRIPLEX_LYS8_18MDA = NeuCodeLys8Triplex.Checked && !mTRAQ.Checked;
            NEUCODE_FOURPLEX_LYS8_12MDA = NeuCodeLys8Fourplex.Checked && !mTRAQ.Checked;
            NEUCODE_SIXPLEX_LYS8_6MDA = NeuCodeLys8Sixplex.Checked && !mTRAQ.Checked;
            NEUCODE_9PLEX_MTRAQ = NeuCodeLys8Triplex.Checked && mTRAQ.Checked;
            NEUCODE_12PLEX_MTRAQ = NeuCodeLys8Fourplex.Checked && mTRAQ.Checked;
            NEUCODE_18PLEX_MTRAQ = NeuCodeLys8Sixplex.Checked && mTRAQ.Checked;
            NEUCODE_DUPLEX_LYS1_6MDA = NeuCodeLys1.Checked;
            NEUCODE_DUPLEX_CARBAMYL = CarbamylCN.Checked;
            NEUCODE_SIXPLEX_MTRAQ = NeuCodeLys8Duplex.Checked && mTRAQ.Checked;
            NEUCODE_SIXPLEX_ARG = NeuCodeLys8Duplex.Checked && Arg.Checked;
            NEUCODE_SIXPLEX_LEU = NeuCodeLys8Duplex.Checked && Leu.Checked;
            SILAC_DUPLEX_LYSC = SILACLys6C.Checked;
            SILAC_DUPLEX_LYSCN = SILACLys8CN.Checked;
            SILAC_DUPLEX_LYSH = SILACLys8D.Checked;
            SILAC_DUPLEX_LEUCN = SILACLeu7CN.Checked;
            SILAC_DUPLEX_LEUH = SILACLeu7D.Checked;
            NEUCODE_12PLEX = Twelveplex.Checked;
            NEUCODE_4PLEX_HEAVY = FourplexH.Checked;
            NEUCODE_4PLEX_MEDIUM = FourplexM.Checked;
            NEUCODE_4PLEX_LIGHT = FourplexL.Checked;
            ICAT = Icat.Checked;

            setExperimentConfiguration();

            WriteMessage("starting");

            // Cycle through .csv files to make a list of identified peptides and properties
            readCsvInputFile(allPeptides, file);   
      
            int rawFileCount = 0;
            int totalRawFiles = RAWFILES.Count;
            int THEORETICALLYRESOLVABLE = 0;
            int CONTAINSLABEL = 0;

            WriteMessage("calculating systematic error");

            // Calculate systematic error by looking for monoisotopes in MS scan preceding best MS/MS scan
            foreach (MSDataFile rawFile in RAWFILES.Values)
            {
                rawFileCount++;
                rawFile.Open();
                WriteMessage("calculating ppm error in raw file " + rawFileCount + " of " + totalRawFiles);
                
                foreach (PeptideID peptide in allPeptides.Values)
                {
                    if (peptide.numLabels > 0)
                    {
                        CONTAINSLABEL++;
                        if (peptide.theoreticallyResolvable)
                        {
                            THEORETICALLYRESOLVABLE++;
                            List<PeptideSpectralMatch> psms = null;
                            if (peptide.PSMs.TryGetValue(rawFile, out psms))
                            {
                                peptide.precursorPPMError(rawFile, PRECURSORPPM);
                            }
                        }
                    }
                }
            }
            
            PRECURSORPPM.Sort();
            //writePrecursorPPMOutput(PRECURSORPPM);
            
            SYSTEMATICERROR = PRECURSORPPM.ElementAt(PRECURSORPPM.Count / 2).Ppm; //Set systematic error as the median value of precursors
            WriteMessage("labeled peptides: " + CONTAINSLABEL);
            WriteMessage("theoretically resolvable peptides: " + THEORETICALLYRESOLVABLE);
            WriteMessage("peptides found: " + (PRECURSORPPM.Count / (double)NUMISOTOPOLOGUES));
            WriteMessage("systematic error: " + SYSTEMATICERROR);
            PRECURSORPPM.Clear();

            WriteMessage("applying systematic error");
            FIRSTSEARCHDONE = true;

            WriteMessage("searching RAW file for pairs");
            rawFileCount = 0;

            foreach (MSDataFile rawFile in RAWFILES.Values)
            {
                rawFileCount++;
                rawFile.Open();
                WriteMessage("finding pairs in raw file " + rawFileCount + " of " + totalRawFiles);

                foreach (PeptideID uniquePeptide in allPeptides.Values)
                {
                    // Only consider peptides that contain at least one label
                    if (uniquePeptide.numLabels > 0 && uniquePeptide.theoreticallyResolvable)
                    {
                        // Only consider peptides that were identified in this raw file
                        List<PeptideSpectralMatch> psms = null;
                        if (uniquePeptide.PSMs.TryGetValue(rawFile, out psms))
                        {
                            uniquePeptide.calculateScanRange(rawFile, RTWINDOW);
                            foreach (MSDataScan currentScan in uniquePeptide.fullScanList)
                            {                               
                                uniquePeptide.findPeaks(currentScan, rawFile);
                            }
                            uniquePeptide.checkPairSpacing(rawFile, spacings);
                        }
                    }
                }
            }

            //writeSpacingOutput(spacings);
            
            if (NOISEBANDCAP && PEAKCOALESCENCE)
            {
                WriteMessage("calculating coalescence threshold");
                MAXIMUMDNL = calculateCoalescenceThreshold(allPeptides, 0.1, INTENSITY_MISSINGCHANNEL);
                WriteMessage("intensity threshold: " + MAXIMUMDNL);
            }

            // Validate NeuCode pairs in which one channel was found to be missing -- apply noise level or discard due to coalescence
            rawFileCount = 0;
            foreach (MSDataFile rawFile in RAWFILES.Values)
            {
                rawFileCount++;
                rawFile.Open();
                WriteMessage("quantifying pairs in raw file " + rawFileCount + " of " + totalRawFiles);

                foreach (PeptideID uniquePeptide in allPeptides.Values)
                {
                    if (uniquePeptide.numLabels > 0 && uniquePeptide.theoreticallyResolvable)
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
                if (uniquePeptide.numLabels > 0 && uniquePeptide.theoreticallyResolvable)
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

            WriteMessage("writing output file");
            writeCsvOutputFile(allPeptides, file);
            WriteMessage("finished");
        }

        public void WriteMessage(string message)
        {
            Console.WriteLine(message);
            richTextBox1.AppendText(message + "\n");
            richTextBox1.ScrollToCaret();
            Application.DoEvents();
        }

        private void writeSpacingOutput(List<Spacing> spacings, string file)
        {
            // Pair spacing outputs can be printed out to a file to assess spacing distributions
            StreamWriter spacingWriter = new StreamWriter(Path.Combine(outputFolderBox.Text, Path.GetFileNameWithoutExtension(file) + "_spacings.csv"));
            string header1 = ("Theo Spacing, Exp Spacing, Charge, MZ");
            spacingWriter.WriteLine(header1);
            WriteMessage("writing output");

            foreach (Spacing spacing in spacings)
            {
                spacingWriter.WriteLine("{0},{1},{2},{3}", spacing.theoSpacing, spacing.spacing, spacing.charge, spacing.MZ);
            }

            spacingWriter.Close();
        }
        
        private void writePrecursorPPMOutput(List<PrecursorPPM> PRECURSORPPM, string file)
        {
            // PPM error outputs can be printed out to a file to assess error distributions
            StreamWriter ppmWriter = new StreamWriter(Path.Combine(outputFolderBox.Text, Path.GetFileNameWithoutExtension(file) + "_ppm.csv"));
            string header = ("Peptide, Charge, PPM Error, E-value");
            ppmWriter.WriteLine(header);
            WriteMessage("writing output");

            foreach (PrecursorPPM ppm in PRECURSORPPM)
            {
                ppmWriter.WriteLine("{0},{1},{2},{3}", ppm.Peptide, ppm.Charge, ppm.Ppm, ppm.EValue);
            }

            ppmWriter.Close();
        }
        
        private void setExperimentConfiguration()
        {
            // 2-plex ICAT
            if (ICAT)
            {
                NEUCODE = false;
                TRADITIONAL = true;
                NUMCHANNELS = 2;
                NUMCLUSTERS = 2;
                NUMISOTOPOLOGUES = 1;
            }

            // 2-plex SILAC
            if (SILAC_DUPLEX_LYSC || SILAC_DUPLEX_LYSCN || SILAC_DUPLEX_LYSH || SILAC_DUPLEX_LEUCN || SILAC_DUPLEX_LEUH)
            {
                NEUCODE = false;
                TRADITIONAL = true;
                NUMCHANNELS = 2;
                NUMCLUSTERS = 2;
                NUMISOTOPOLOGUES = 1;

                if (IncompleteIncorporation.Checked)
                {
                    CHECKPARTIALINCORPORATION = true;
                }

                if (SILAC_DUPLEX_LEUCN && Conversion.Checked)
                {
                    CORRECTLEUNLOSS = true;
                }
            }

            // 2-plex NeuCode
            if (NEUCODE_DUPLEX_LYS8_36MDA || NEUCODE_DUPLEX_LYS1_6MDA || NEUCODE_DUPLEX_CARBAMYL)
            {
                NEUCODE = true;
                TRADITIONAL = false;
                NUMCHANNELS = 2;
                NUMCLUSTERS = 1;
                NUMISOTOPOLOGUES = 2;
            }

            // 3-plex NeuCode
            if (NEUCODE_TRIPLEX_LYS8_18MDA)
            {
                NEUCODE = true;
                TRADITIONAL = false;
                NUMCHANNELS = 3;
                NUMCLUSTERS = 1;
                NUMISOTOPOLOGUES = 3;
            }

            // 4-plex NeuCode
            if (NEUCODE_FOURPLEX_LYS8_12MDA || NEUCODE_4PLEX_LIGHT || NEUCODE_4PLEX_MEDIUM || NEUCODE_4PLEX_HEAVY)
            {
                NEUCODE = true;
                TRADITIONAL = false;
                NUMCHANNELS = 4;
                NUMCLUSTERS = 1;
                NUMISOTOPOLOGUES = 4;
            }

            // 6-plex NeuCode
            if (NEUCODE_SIXPLEX_LYS8_6MDA || NEUCODE_SIXPLEX_MTRAQ || NEUCODE_SIXPLEX_ARG || NEUCODE_SIXPLEX_LEU)
            {
                NEUCODE = true;
                TRADITIONAL = false;
                NUMCHANNELS = 6;

                if (NEUCODE_SIXPLEX_MTRAQ || NEUCODE_SIXPLEX_ARG || NEUCODE_SIXPLEX_LEU)
                {
                    NUMCLUSTERS = 3;
                    NUMISOTOPOLOGUES = 2;

                    if (NEUCODE_SIXPLEX_ARG && Conversion.Checked)
                    {
                        CORRECTARGPROLINE = true;
                    }

                    if (NEUCODE_SIXPLEX_LEU && Conversion.Checked)
                    {
                        CORRECTLEUDLOSS = true;
                    }
                }
                else
                {
                    NUMCLUSTERS = 1;
                    NUMISOTOPOLOGUES = 6;

                    if (AGCBINNING)
                    {
                        SEGMENTEDINJECTIONTIMES = true;
                        INJECTIONTIMES = new Dictionary<string, Dictionary<int, Dictionary<Range<double>, double>>>();
                    }
                }
            }

            // 9-plex NeuCode
            if (NEUCODE_9PLEX_MTRAQ)
            {
                NEUCODE = true;
                TRADITIONAL = false;
                NUMCHANNELS = 9;
                NUMCLUSTERS = 3;
                NUMISOTOPOLOGUES = 3;
            }

            // 12-plex NeuCode
            if (NEUCODE_12PLEX || NEUCODE_12PLEX_MTRAQ)
            {
                NEUCODE = true;
                TRADITIONAL = false;
                NUMCHANNELS = 12;
                NUMCLUSTERS = 3;
                NUMISOTOPOLOGUES = 4;
            }

            // 18-plex NeuCode
            if (NEUCODE_18PLEX_MTRAQ)
            {
                NEUCODE = true;
                TRADITIONAL = false;
                NUMCHANNELS = 18;
                NUMCLUSTERS = 3;
                NUMISOTOPOLOGUES = 6;
            }
        }
        
        private void writeCsvOutputFile(Dictionary<string, PeptideID> allPeptides, string file)
        {
            string outputName;
            if (NEUCODE)
            {
                outputName = Path.Combine(outputFolderBox.Text, Path.GetFileNameWithoutExtension(file) + "_NeuCode_Quant.csv");
            }
            else
            {
                outputName = Path.Combine(outputFolderBox.Text, Path.GetFileNameWithoutExtension(file) + "_SILAC_Quant.csv");
            }
            StreamWriter writer1 = new StreamWriter(outputName);

            if (NUMCHANNELS == 2)
            {
                string header1;
                if (SILAC_DUPLEX_LEUCN || SILAC_DUPLEX_LEUH)
                {
                    header1 = ("Scan number(s), Raw File(s), Charge State(s), Best Scan Number, Peptide, # Labels, Theo Mass 1, Theo Mass 2, Adjusted Mass 1, Adjusted Mass 2, Label Conversion?, Quantified no NBC MS1 Scans, # no NBC Measurements, No NBC Intensity 1, No NBC Intensity 2, Count 1, Count 2");
                    writer1.WriteLine(header1);
                    foreach (PeptideID peptide in allPeptides.Values)
                    {
                        writer1.WriteLine("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16}",
                            peptide.scanNumbers, peptide.rawFiles, peptide.chargeStates, peptide.bestPSM.ScanNumber, peptide.sequence, peptide.numLabels,
                            peptide.theoMasses[0, 0], peptide.theoMasses[1, 0],
                            peptide.adjustedTheoMasses[0, 0], peptide.adjustedTheoMasses[1, 0], peptide.conversionFactor, 
                            peptide.countCompletePairs, peptide.countCompleteIsotopes[0], peptide.completeTotalIntensity[0, NUMISOTOPES], peptide.completeTotalIntensity[1, NUMISOTOPES], peptide.finalQuantified[0], peptide.finalQuantified[1]);
                    }
                }
                else if (!NOISEBANDCAP)
                {
                    header1 = ("Scan number(s), Raw File(s), Charge State(s), Best Scan Number, Peptide, # Labels, Theo Mass 1, Theo Mass 2, Adjusted Mass 1, Adjusted Mass 2, Resolvable?, Quantified MS1 Scans, # Measurements, Intensity 1, Intensity 2, RT 1, RT 2, Ratio 2/1, Ratio Count, Missing Channels?");
                    writer1.WriteLine(header1);
                    foreach (PeptideID peptide in allPeptides.Values)
                    {
                        writer1.WriteLine("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19}",
                            peptide.scanNumbers, peptide.rawFiles, peptide.chargeStates, peptide.bestPSM.ScanNumber, peptide.sequence, peptide.numLabels,
                            peptide.theoMasses[0, 0], peptide.theoMasses[1, 0],
                            peptide.adjustedTheoMasses[0, 0], peptide.adjustedTheoMasses[1, 0], peptide.GetTheoreticalResolvability(0),
                            peptide.countCompletePairs, peptide.countCompleteIsotopes[0], peptide.completeTotalIntensity[0, NUMISOTOPES], peptide.completeTotalIntensity[1, NUMISOTOPES], 
                            peptide.maxCompleteIntensity[0,1], peptide.maxCompleteIntensity[1,1],
                            peptide.heavyToLightRatioSum[0, 0], peptide.finalQuantified[0], peptide.quantifiedNoiseIncluded[0]);
                    }
                }
                else
                {
                    header1 = ("Scan number(s), Raw File(s), Charge State(s), Best Scan Number, Peptide, # Labels, Theo Mass 1, Theo Mass 2, Adjusted Mass 1, Adjusted Mass 2, Resolvable?, Quantified MS1 Scans, # Measurements, Intensity 1, Intensity 2, Ratio 2/1, Ratio Count, Missing Channels?");
                    writer1.WriteLine(header1);
                    foreach (PeptideID peptide in allPeptides.Values)
                    {
                        writer1.WriteLine("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17}",
                            peptide.scanNumbers, peptide.rawFiles, peptide.chargeStates, peptide.bestPSM.ScanNumber, peptide.sequence, peptide.numLabels,
                            peptide.theoMasses[0, 0], peptide.theoMasses[1, 0],
                            peptide.adjustedTheoMasses[0, 0], peptide.adjustedTheoMasses[1, 0], peptide.GetTheoreticalResolvability(0),
                            peptide.countAllPairs, peptide.countAllIsotopes[0], peptide.totalIntensity[0, NUMISOTOPES], peptide.totalIntensity[1, NUMISOTOPES], 
                            peptide.heavyToLightRatioSum[0, 0], peptide.finalQuantified[0], peptide.quantifiedNoiseIncluded[0]);
                    }
                }
                writer1.Close();
            }
            else if (NUMCHANNELS == 3)
            {
                string header1;
                if (!NOISEBANDCAP)
                {
                    header1 = ("Scan number(s), Raw File(s), Charge State(s), Best Scan Number, Peptide, # Labels, Theo Mass 1, Theo Mass 2, Theo Mass 3, Adjusted Mass 1, Adjusted Mass 2, Adjusted Mass 3, Resolvable?, Quantified MS1 Scans, # Measurements, Intensity 1, Intensity 2, Intensity 3, RT 1, RT 2, RT 3, Ratio 2/1, Ratio 3/1, Ratio Count, Missing Channels?");
                    writer1.WriteLine(header1);
                    foreach (PeptideID peptide in allPeptides.Values)
                    {
                        writer1.WriteLine("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23},{24}",
                            peptide.scanNumbers, peptide.rawFiles, peptide.chargeStates, peptide.bestPSM.ScanNumber, peptide.sequence, peptide.numLabels,
                            peptide.theoMasses[0, 0], peptide.theoMasses[1, 0], peptide.theoMasses[2,0],
                            peptide.adjustedTheoMasses[0, 0], peptide.adjustedTheoMasses[1, 0], peptide.adjustedTheoMasses[2,0], peptide.GetTheoreticalResolvability(0),
                            peptide.countCompletePairs, peptide.countCompleteIsotopes[0], peptide.completeTotalIntensity[0, NUMISOTOPES], peptide.completeTotalIntensity[1, NUMISOTOPES], peptide.completeTotalIntensity[2,NUMISOTOPES],
                            peptide.maxCompleteIntensity[0,1], peptide.maxCompleteIntensity[1,1], peptide.maxCompleteIntensity[2,1],
                            peptide.heavyToLightRatioSum[0, 0], peptide.heavyToLightRatioSum[1,0], peptide.finalQuantified[0], peptide.quantifiedNoiseIncluded[0]);
                    }
                }
                else
                {
                    header1 = ("Scan number(s), Raw File(s), Charge State(s), Best Scan Number, Peptide, # Labels, Theo Mass 1, Theo Mass 2, Theo Mass 3, Adjusted Mass 1, Adjusted Mass 2, Adjusted Mass 3, Resolvable?, Quantified MS1 Scans, # Measurements, Intensity 1, Intensity 2, Intensity 3, Ratio 2/1, Ratio 3/1, Ratio Count, Missing Channels?");
                    writer1.WriteLine(header1);
                    foreach (PeptideID peptide in allPeptides.Values)
                    {
                        writer1.WriteLine("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21}",
                            peptide.scanNumbers, peptide.rawFiles, peptide.chargeStates, peptide.bestPSM.ScanNumber, peptide.sequence, peptide.numLabels,
                            peptide.theoMasses[0, 0], peptide.theoMasses[1, 0], peptide.theoMasses[2, 0],
                            peptide.adjustedTheoMasses[0, 0], peptide.adjustedTheoMasses[1, 0], peptide.adjustedTheoMasses[2, 0], peptide.GetTheoreticalResolvability(0),
                            peptide.countAllPairs, peptide.countAllIsotopes[0], peptide.totalIntensity[0, NUMISOTOPES], peptide.totalIntensity[1, NUMISOTOPES], peptide.totalIntensity[2, NUMISOTOPES], 
                            peptide.heavyToLightRatioSum[0, 0], peptide.heavyToLightRatioSum[1, 0], peptide.finalQuantified[0], peptide.quantifiedNoiseIncluded[0]);
                    }
                }
                writer1.Close();
            }
            else if (NUMCHANNELS == 4)
            {
                string header1;
                if (!NOISEBANDCAP)
                {
                    header1 = ("Scan number(s), Raw File(s), Charge State(s), Best Scan Number, Peptide, # Labels, Theo Mass 1, Theo Mass 2, Theo Mass 3, Theo Mass 4, Adjusted Mass 1, Adjusted Mass 2, Adjusted Mass 3, Adjusted Mass 4, Resolvable?, Quantified MS1 Scans, # Measurements, Intensity 1, Intensity 2, Intensity 3, Intensity 4, RT 1, RT 2, RT 3, RT 4, Ratio 2/1, Ratio 3/1, Ratio 4/1, Ratio Count, Missing Channels?");
                    writer1.WriteLine(header1);
                    foreach (PeptideID peptide in allPeptides.Values)
                    {
                        writer1.WriteLine("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23},{24},{25},{26},{27},{28},{29}",
                            peptide.scanNumbers, peptide.rawFiles, peptide.chargeStates, peptide.bestPSM.ScanNumber, peptide.sequence, peptide.numLabels,
                            peptide.theoMasses[0, 0], peptide.theoMasses[1, 0], peptide.theoMasses[2, 0], peptide.theoMasses[3,0],
                            peptide.adjustedTheoMasses[0, 0], peptide.adjustedTheoMasses[1, 0], peptide.adjustedTheoMasses[2, 0], peptide.adjustedTheoMasses[3,0], peptide.GetTheoreticalResolvability(0),
                            peptide.countCompletePairs, peptide.countCompleteIsotopes[0], peptide.completeTotalIntensity[0, NUMISOTOPES], peptide.completeTotalIntensity[1, NUMISOTOPES], peptide.completeTotalIntensity[2, NUMISOTOPES], peptide.completeTotalIntensity[3,NUMISOTOPES], 
                            peptide.maxCompleteIntensity[0,1], peptide.maxCompleteIntensity[1,1], peptide.maxCompleteIntensity[2,1], peptide.maxCompleteIntensity[3,1],
                            peptide.heavyToLightRatioSum[0, 0], peptide.heavyToLightRatioSum[1, 0], peptide.heavyToLightRatioSum[2,0], peptide.finalQuantified[0], peptide.quantifiedNoiseIncluded[0]);
                    }
                }
                else
                {
                    header1 = ("Scan number(s), Raw File(s), Charge State(s), Best Scan Number, Peptide, # Labels, Theo Mass 1, Theo Mass 2, Theo Mass 3, Theo Mass 4, Adjusted Mass 1, Adjusted Mass 2, Adjusted Mass 3, Adjusted Mass 4, Resolvable?, Quantified MS1 Scans, # Measurements, Intensity 1, Intensity 2, Intensity 3, Intensity 4, Ratio 2/1, Ratio 3/1, Ratio 4/1, Ratio Count, Missing Channels?");
                    writer1.WriteLine(header1);
                    foreach (PeptideID peptide in allPeptides.Values)
                    {
                        writer1.WriteLine("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23},{24},{25}",
                            peptide.scanNumbers, peptide.rawFiles, peptide.chargeStates, peptide.bestPSM.ScanNumber, peptide.sequence, peptide.numLabels,
                            peptide.theoMasses[0, 0], peptide.theoMasses[1, 0], peptide.theoMasses[2, 0], peptide.theoMasses[3, 0],
                            peptide.adjustedTheoMasses[0, 0], peptide.adjustedTheoMasses[1, 0], peptide.adjustedTheoMasses[2, 0], peptide.adjustedTheoMasses[3, 0], peptide.GetTheoreticalResolvability(0),
                            peptide.countAllPairs, peptide.countAllIsotopes[0], peptide.totalIntensity[0, NUMISOTOPES], peptide.totalIntensity[1, NUMISOTOPES], peptide.totalIntensity[2, NUMISOTOPES], peptide.totalIntensity[3, NUMISOTOPES], 
                            peptide.heavyToLightRatioSum[0, 0], peptide.heavyToLightRatioSum[1, 0], peptide.heavyToLightRatioSum[2, 0], peptide.finalQuantified[0], peptide.quantifiedNoiseIncluded[0]);
                    }
                }
                writer1.Close();
            }
            else if (NUMCHANNELS == 6)
            {
                string header1;
                if (!NOISEBANDCAP)
                {
                    if (NEUCODE_SIXPLEX_LYS8_6MDA)
                    {
                        header1 = ("Scan number(s), Raw File(s), Charge State(s), Best Scan Number, Peptide, # Labels, Theo Mass 1, Theo Mass 2, Theo Mass 3, Theo Mass 4, Theo Mass 5, Theo Mass 6, Adjusted Mass 1, Adjusted Mass 2, Adjusted Mass 3, Adjusted Mass 4, Adjusted Mass 5, Adjusted Mass 6, Resolvable?, Quantified MS1 Scans, # Measurements, Intensity 1, Intensity 2, Intensity 3, Intensity 4, Intensity 5, Intensity 6, RT 1, RT 2, RT 3, RT 4, RT 5, RT 6, Ratio 2/1, Ratio 3/1, Ratio 4/1, Ratio 5/1, Ratio 6/1, Ratio Count, Missing Channels?");
                        writer1.WriteLine(header1);
                        foreach (PeptideID peptide in allPeptides.Values)
                        {
                            writer1.WriteLine("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23},{24},{25},{26},{27},{28},{29},{30},{31},{32},{33},{34},{35},{36},{37},{38},{39}",
                                peptide.scanNumbers, peptide.rawFiles, peptide.chargeStates, peptide.bestPSM.ScanNumber, peptide.sequence, peptide.numLabels,
                                peptide.theoMasses[0, 0], peptide.theoMasses[1, 0], peptide.theoMasses[2, 0], peptide.theoMasses[3, 0], peptide.theoMasses[4,0], peptide.theoMasses[5,0],
                                peptide.adjustedTheoMasses[0, 0], peptide.adjustedTheoMasses[1, 0], peptide.adjustedTheoMasses[2, 0], peptide.adjustedTheoMasses[3, 0], peptide.adjustedTheoMasses[4,0], peptide.adjustedTheoMasses[5,0], peptide.GetTheoreticalResolvability(0),
                                peptide.countCompletePairs, peptide.countCompleteIsotopes[0], peptide.completeTotalIntensity[0, NUMISOTOPES], peptide.completeTotalIntensity[1, NUMISOTOPES], peptide.completeTotalIntensity[2, NUMISOTOPES], peptide.completeTotalIntensity[3, NUMISOTOPES], peptide.completeTotalIntensity[4, NUMISOTOPES], peptide.completeTotalIntensity[5, NUMISOTOPES], 
                                peptide.maxCompleteIntensity[0,1], peptide.maxCompleteIntensity[1,1], peptide.maxCompleteIntensity[2,1], peptide.maxCompleteIntensity[3,1], peptide.maxCompleteIntensity[4,1], peptide.maxCompleteIntensity[5,1],
                                peptide.heavyToLightRatioSum[0, 0], peptide.heavyToLightRatioSum[1, 0], peptide.heavyToLightRatioSum[2, 0], peptide.heavyToLightRatioSum[3,0], peptide.heavyToLightRatioSum[4,0], peptide.finalQuantified[0], peptide.quantifiedNoiseIncluded[0]);
                        }
                    }
                    else if (NEUCODE_SIXPLEX_MTRAQ || NEUCODE_SIXPLEX_ARG)
                    {
                        header1 = ("Scan number(s), Raw File(s), Charge State(s), Best Scan Number, Peptide, # Labels, Theo Mass 1, Theo Mass 2, Theo Mass 3, Theo Mass 4, Theo Mass 5, Theo Mass 6, Adjusted Mass 1, Adjusted Mass 2, Adjusted Mass 3, Adjusted Mass 4, Adjusted Mass 5, Adjusted Mass 6, Cluster 1 Resolvable?, Cluster 2 Resolvable?, Cluster 3 Resolvable?, Quantified MS1 Scans, Cluster 1 # Measurements, Cluster 2 # Measurements, Cluster 3 # Measurements, Intensity 1, Intensity 2, Intensity 3, Intensity 4, Intensity 5, Intensity 6, Ratio 2/1, Ratio 3/4, Ratio 5/6, Cluster 1 Ratio Count, Cluster 2 Ratio Count, Cluster 3 Ratio Count, Cluster 1 Missing Channels?, Cluster 2 Missing Channels?, Cluster 3 Missing Channels?");
                        writer1.WriteLine(header1);
                        foreach (PeptideID peptide in allPeptides.Values)
                        {
                            writer1.WriteLine("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23},{24},{25},{26},{27},{28},{29},{30},{31},{32},{33},{34},{35},{36},{37},{38},{39}",
                                peptide.scanNumbers, peptide.rawFiles, peptide.chargeStates, peptide.bestPSM.ScanNumber, peptide.sequence, peptide.numLabels,
                                peptide.theoMasses[0, 0], peptide.theoMasses[1, 0], peptide.theoMasses[2, 0], peptide.theoMasses[3, 0], peptide.theoMasses[4, 0], peptide.theoMasses[5, 0],
                                peptide.adjustedTheoMasses[0, 0], peptide.adjustedTheoMasses[1, 0], peptide.adjustedTheoMasses[2, 0], peptide.adjustedTheoMasses[3, 0], peptide.adjustedTheoMasses[4, 0], peptide.adjustedTheoMasses[5, 0], peptide.GetTheoreticalResolvability(0), peptide.GetTheoreticalResolvability(1), peptide.GetTheoreticalResolvability(2),
                                peptide.countCompletePairs, peptide.countCompleteIsotopes[0], peptide.countCompleteIsotopes[1], peptide.countCompleteIsotopes[2], peptide.completeTotalIntensity[0, NUMISOTOPES], peptide.completeTotalIntensity[1, NUMISOTOPES], peptide.completeTotalIntensity[2, NUMISOTOPES], peptide.completeTotalIntensity[3, NUMISOTOPES], peptide.completeTotalIntensity[4, NUMISOTOPES], peptide.completeTotalIntensity[5, NUMISOTOPES], 
                                peptide.heavyToLightRatioSum[0, 0], peptide.heavyToLightRatioSum[1, 0], peptide.heavyToLightRatioSum[2, 0], peptide.finalQuantified[0], peptide.finalQuantified[1], peptide.finalQuantified[2], peptide.quantifiedNoiseIncluded[0], peptide.quantifiedNoiseIncluded[1], peptide.quantifiedNoiseIncluded[2]);
                        }
                    }
                }
                else
                {
                    if (NEUCODE_SIXPLEX_LYS8_6MDA)
                    {
                        header1 = ("Scan number(s), Raw File(s), Charge State(s), Best Scan Number, Peptide, # Labels, Theo Mass 1, Theo Mass 2, Theo Mass 3, Theo Mass 4, Theo Mass 5, Theo Mass 6, Adjusted Mass 1, Adjusted Mass 2, Adjusted Mass 3, Adjusted Mass 4, Adjusted Mass 5, Adjusted Mass 6, Resolvable?, Quantified MS1 Scans, # Measurements, Intensity 1, Intensity 2, Intensity 3, Intensity 4, Intensity 5, Intensity 6, Ratio 2/1, Ratio 3/1, Ratio 4/1, Ratio 5/1, Ratio 6/1, Ratio Count, Missing Channels?");
                        writer1.WriteLine(header1);
                        foreach (PeptideID peptide in allPeptides.Values)
                        {
                            writer1.WriteLine("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23},{24},{25},{26},{27},{28},{29},{30},{31},{32},{33}",
                                peptide.scanNumbers, peptide.rawFiles, peptide.chargeStates, peptide.bestPSM.ScanNumber, peptide.sequence, peptide.numLabels,
                                peptide.theoMasses[0, 0], peptide.theoMasses[1, 0], peptide.theoMasses[2, 0], peptide.theoMasses[3, 0], peptide.theoMasses[4, 0], peptide.theoMasses[5, 0],
                                peptide.adjustedTheoMasses[0, 0], peptide.adjustedTheoMasses[1, 0], peptide.adjustedTheoMasses[2, 0], peptide.adjustedTheoMasses[3, 0], peptide.adjustedTheoMasses[4, 0], peptide.adjustedTheoMasses[5, 0], peptide.GetTheoreticalResolvability(0),
                                peptide.countAllPairs, peptide.countAllIsotopes[0], peptide.totalIntensity[0, NUMISOTOPES], peptide.totalIntensity[1, NUMISOTOPES], peptide.totalIntensity[2, NUMISOTOPES], peptide.totalIntensity[3, NUMISOTOPES], peptide.totalIntensity[4, NUMISOTOPES], peptide.totalIntensity[5, NUMISOTOPES],
                                peptide.heavyToLightRatioSum[0, 0], peptide.heavyToLightRatioSum[1, 0], peptide.heavyToLightRatioSum[2, 0], peptide.heavyToLightRatioSum[3, 0], peptide.heavyToLightRatioSum[4, 0], peptide.finalQuantified[0], peptide.quantifiedNoiseIncluded[0]);
                        }
                    }
                    else if (NEUCODE_SIXPLEX_MTRAQ || NEUCODE_SIXPLEX_ARG)
                    {
                        header1 = ("Scan number(s), Raw File(s), Charge State(s), Best Scan Number, Peptide, # Labels, Theo Mass 1, Theo Mass 2, Theo Mass 3, Theo Mass 4, Theo Mass 5, Theo Mass 6, Adjusted Mass 1, Adjusted Mass 2, Adjusted Mass 3, Adjusted Mass 4, Adjusted Mass 5, Adjusted Mass 6, Cluster 1 Resolvable?, Cluster 2 Resolvable?, Cluster 3 Resolvable?, Quantified MS1 Scans, Cluster 1 # Measurements, Cluster 2 # Measurements, Cluster 3 # Measurements, Intensity 1, Intensity 2, Intensity 3, Intensity 4, Intensity 5, Intensity 6, Ratio 2/1, Ratio 3/4, Ratio 5/6, Cluster 1 Ratio Count, Cluster 2 Ratio Count, Cluster 3 Ratio Count, Cluster 1 Missing Channels?, Cluster 2 Missing Channels?, Cluster 3 Missing Channels?");
                        writer1.WriteLine(header1);
                        foreach (PeptideID peptide in allPeptides.Values)
                        {
                            writer1.WriteLine("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23},{24},{25},{26},{27},{28},{29},{30},{31},{32},{33},{34},{35},{36},{37},{38},{39}",
                                peptide.scanNumbers, peptide.rawFiles, peptide.chargeStates, peptide.bestPSM.ScanNumber, peptide.sequence, peptide.numLabels,
                                peptide.theoMasses[0, 0], peptide.theoMasses[1, 0], peptide.theoMasses[2, 0], peptide.theoMasses[3, 0], peptide.theoMasses[4, 0], peptide.theoMasses[5, 0],
                                peptide.adjustedTheoMasses[0, 0], peptide.adjustedTheoMasses[1, 0], peptide.adjustedTheoMasses[2, 0], peptide.adjustedTheoMasses[3, 0], peptide.adjustedTheoMasses[4, 0], peptide.adjustedTheoMasses[5, 0], peptide.GetTheoreticalResolvability(0), peptide.GetTheoreticalResolvability(1), peptide.GetTheoreticalResolvability(2),
                                peptide.countAllPairs, peptide.countAllIsotopes[0], peptide.countAllIsotopes[1], peptide.countAllIsotopes[2], peptide.totalIntensity[0, NUMISOTOPES], peptide.totalIntensity[1, NUMISOTOPES], peptide.totalIntensity[2, NUMISOTOPES], peptide.totalIntensity[3, NUMISOTOPES], peptide.totalIntensity[4, NUMISOTOPES], peptide.totalIntensity[5, NUMISOTOPES],
                                peptide.heavyToLightRatioSum[0, 0], peptide.heavyToLightRatioSum[1, 0], peptide.heavyToLightRatioSum[2, 0], peptide.finalQuantified[0], peptide.finalQuantified[1], peptide.finalQuantified[2], peptide.quantifiedNoiseIncluded[0], peptide.quantifiedNoiseIncluded[1], peptide.quantifiedNoiseIncluded[2]);
                        }
                    }
                }
                writer1.Close();
            }
            else if (NUMCHANNELS == 9)
            {
                string header1;
                if (!NOISEBANDCAP)
                {
                    header1 = ("Scan number(s), Raw File(s), Charge State(s), Best Scan Number, Peptide, # Labels, Theo Mass 1, Theo Mass 2, Theo Mass 3, Theo Mass 4, Theo Mass 5, Theo Mass 6, Theo Mass 7, Theo Mass 8, Theo Mass 9, Adjusted Mass 1, Adjusted Mass 2, Adjusted Mass 3, Adjusted Mass 4, Adjusted Mass 5, Adjusted Mass 6, Adjusted Mass 7, Adjusted Mass 8, Adjusted Mass 9, Cluster 1 Resolvable?, Cluster 2 Resolvable?, Cluster 3 Resolvable?, Quantified MS1 Scans, Cluster 1 # Measurements, Cluster 2 # Measurements, Cluster 3 # Measurements, Intensity 1, Intensity 2, Intensity 3, Intensity 4, Intensity 5, Intensity 6, Intensity 7, Intensity 8, Intensity 9, Ratio 2/1, Ratio 3/1, Ratio 5/4, Ratio 6/4, Ratio 8/7, Ratio 9/7, Cluster 1 Ratio Count, Cluster 2 Ratio Count, Cluster 3 Ratio Count, Cluster 1 Missing Channels?, Cluster 2 Missing Channels?, Cluster 3 Missing Channels?");
                    writer1.WriteLine(header1);
                    foreach (PeptideID peptide in allPeptides.Values)
                    {
                        writer1.WriteLine("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23},{24},{25},{26},{27},{28},{29},{30},{31},{32},{33},{34},{35},{36},{37},{38},{39},{40},{41},{42},{43},{44},{45},{46},{47},{48},{49},{50},{51}",
                            peptide.scanNumbers, peptide.rawFiles, peptide.chargeStates, peptide.bestPSM.ScanNumber, peptide.sequence, peptide.numLabels,
                            peptide.theoMasses[0, 0], peptide.theoMasses[1, 0], peptide.theoMasses[2, 0], peptide.theoMasses[3, 0], peptide.theoMasses[4, 0], peptide.theoMasses[5, 0], peptide.theoMasses[6,0], peptide.theoMasses[7,0], peptide.theoMasses[8,0],
                            peptide.adjustedTheoMasses[0, 0], peptide.adjustedTheoMasses[1, 0], peptide.adjustedTheoMasses[2, 0], peptide.adjustedTheoMasses[3, 0], peptide.adjustedTheoMasses[4, 0], peptide.adjustedTheoMasses[5, 0], peptide.adjustedTheoMasses[6,0], peptide.adjustedTheoMasses[7,0], peptide.adjustedTheoMasses[8,0], peptide.GetTheoreticalResolvability(0), peptide.GetTheoreticalResolvability(1), peptide.GetTheoreticalResolvability(2),
                            peptide.countCompletePairs, peptide.countCompleteIsotopes[0], peptide.countCompleteIsotopes[1], peptide.countCompleteIsotopes[2], peptide.completeTotalIntensity[0, NUMISOTOPES], peptide.completeTotalIntensity[1, NUMISOTOPES], peptide.completeTotalIntensity[2, NUMISOTOPES], peptide.completeTotalIntensity[3, NUMISOTOPES], peptide.completeTotalIntensity[4, NUMISOTOPES], peptide.completeTotalIntensity[5, NUMISOTOPES], peptide.completeTotalIntensity[6, NUMISOTOPES], peptide.completeTotalIntensity[7, NUMISOTOPES], peptide.completeTotalIntensity[8, NUMISOTOPES],
                            peptide.heavyToLightRatioSum[0, 0], peptide.heavyToLightRatioSum[1, 0], peptide.heavyToLightRatioSum[2, 0], peptide.heavyToLightRatioSum[3,0], peptide.heavyToLightRatioSum[4,0], peptide.heavyToLightRatioSum[5,0], peptide.finalQuantified[0], peptide.finalQuantified[1], peptide.finalQuantified[2], peptide.quantifiedNoiseIncluded[0], peptide.quantifiedNoiseIncluded[1], peptide.quantifiedNoiseIncluded[2]);
                    }
                }
                else
                {
                    header1 = ("Scan number(s), Raw File(s), Charge State(s), Best Scan Number, Peptide, # Labels, Theo Mass 1, Theo Mass 2, Theo Mass 3, Theo Mass 4, Theo Mass 5, Theo Mass 6, Theo Mass 7, Theo Mass 8, Theo Mass 9, Adjusted Mass 1, Adjusted Mass 2, Adjusted Mass 3, Adjusted Mass 4, Adjusted Mass 5, Adjusted Mass 6, Adjusted Mass 7, Adjusted Mass 8, Adjusted Mass 9, Cluster 1 Resolvable?, Cluster 2 Resolvable?, Cluster 3 Resolvable?, Quantified MS1 Scans, Cluster 1 # Measurements, Cluster 2 # Measurements, Cluster 3 # Measurements, Intensity 1, Intensity 2, Intensity 3, Intensity 4, Intensity 5, Intensity 6, Intensity 7, Intensity 8, Intensity 9, Ratio 2/1, Ratio 3/1, Ratio 5/4, Ratio 6/4, Ratio 8/7, Ratio 9/7, Cluster 1 Ratio Count, Cluster 2 Ratio Count, Cluster 3 Ratio Count, Cluster 1 Missing Channels?, Cluster 2 Missing Channels?, Cluster 3 Missing Channels?");
                    writer1.WriteLine(header1);
                    foreach (PeptideID peptide in allPeptides.Values)
                    {
                        writer1.WriteLine("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23},{24},{25},{26},{27},{28},{29},{30},{31},{32},{33},{34},{35},{36},{37},{38},{39},{40},{41},{42},{43},{44},{45},{46},{47},{48},{49},{50},{51}",
                            peptide.scanNumbers, peptide.rawFiles, peptide.chargeStates, peptide.bestPSM.ScanNumber, peptide.sequence, peptide.numLabels,
                            peptide.theoMasses[0, 0], peptide.theoMasses[1, 0], peptide.theoMasses[2, 0], peptide.theoMasses[3, 0], peptide.theoMasses[4, 0], peptide.theoMasses[5, 0], peptide.theoMasses[6, 0], peptide.theoMasses[7, 0], peptide.theoMasses[8, 0],
                            peptide.adjustedTheoMasses[0, 0], peptide.adjustedTheoMasses[1, 0], peptide.adjustedTheoMasses[2, 0], peptide.adjustedTheoMasses[3, 0], peptide.adjustedTheoMasses[4, 0], peptide.adjustedTheoMasses[5, 0], peptide.adjustedTheoMasses[6, 0], peptide.adjustedTheoMasses[7, 0], peptide.adjustedTheoMasses[8, 0], peptide.GetTheoreticalResolvability(0), peptide.GetTheoreticalResolvability(1), peptide.GetTheoreticalResolvability(2),
                            peptide.countAllPairs, peptide.countAllIsotopes[0], peptide.countAllIsotopes[1], peptide.countAllIsotopes[2], peptide.totalIntensity[0, NUMISOTOPES], peptide.totalIntensity[1, NUMISOTOPES], peptide.totalIntensity[2, NUMISOTOPES], peptide.totalIntensity[3, NUMISOTOPES], peptide.totalIntensity[4, NUMISOTOPES], peptide.totalIntensity[5, NUMISOTOPES], peptide.totalIntensity[6, NUMISOTOPES], peptide.totalIntensity[7, NUMISOTOPES], peptide.totalIntensity[8, NUMISOTOPES],
                            peptide.heavyToLightRatioSum[0, 0], peptide.heavyToLightRatioSum[1, 0], peptide.heavyToLightRatioSum[2, 0], peptide.heavyToLightRatioSum[3, 0], peptide.heavyToLightRatioSum[4, 0], peptide.heavyToLightRatioSum[5, 0], peptide.finalQuantified[0], peptide.finalQuantified[1], peptide.finalQuantified[2], peptide.quantifiedNoiseIncluded[0], peptide.quantifiedNoiseIncluded[1], peptide.quantifiedNoiseIncluded[2]);
                    }
                }
                writer1.Close();
            }
            else if (NUMCHANNELS == 12)
            {
                string header1;
                if (!NOISEBANDCAP)
                {
                    header1 = ("Scan number(s), Raw File(s), Charge State(s), Best Scan Number, Peptide, # Labels, Theo Mass 1, Theo Mass 2, Theo Mass 3, Theo Mass 4, Theo Mass 5, Theo Mass 6, Theo Mass 7, Theo Mass 8, Theo Mass 9, Theo Mass 10, Theo Mass 11, Theo Mass 12, Adjusted Mass 1, Adjusted Mass 2, Adjusted Mass 3, Adjusted Mass 4, Adjusted Mass 5, Adjusted Mass 6, Adjusted Mass 7, Adjusted Mass 8, Adjusted Mass 9, Adjusted Mass 10, Adjusted Mass 11, Adjusted Mass 12, Cluster 1 Resolvable?, Cluster 2 Resolvable?, Cluster 3 Resolvable?, Quantified MS1 Scans, Cluster 1 # Measurements, Cluster 2 # Measurements, Cluster 3 # Measurements, Intensity 1, Intensity 2, Intensity 3, Intensity 4, Intensity 5, Intensity 6, Intensity 7, Intensity 8, Intensity 9, Intensity 10, Intensity 11, Intensity 12, Ratio 2/1, Ratio 3/1, Ratio 4/1, Ratio 6/5, Ratio 7/5, Ratio 8/5, Ratio 10/9, Ratio 11/9, Ratio 12/9, Cluster 1 Ratio Count, Cluster 2 Ratio Count, Cluster 3 Ratio Count, Cluster 1 Missing Channels?, Cluster 2 Missing Channels?, Cluster 3 Missing Channels?");
                    writer1.WriteLine(header1);
                    foreach (PeptideID peptide in allPeptides.Values)
                    {
                        writer1.WriteLine("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23},{24},{25},{26},{27},{28},{29},{30},{31},{32},{33},{34},{35},{36},{37},{38},{39},{40},{41},{42},{43},{44},{45},{46},{47},{48},{49},{50},{51},{52},{53},{54},{55},{56},{57},{58},{59},{60},{61},{62},{63}",
                            peptide.scanNumbers, peptide.rawFiles, peptide.chargeStates, peptide.bestPSM.ScanNumber, peptide.sequence, peptide.numLabels,
                            peptide.theoMasses[0, 0], peptide.theoMasses[1, 0], peptide.theoMasses[2, 0], peptide.theoMasses[3, 0], peptide.theoMasses[4, 0], peptide.theoMasses[5, 0], peptide.theoMasses[6, 0], peptide.theoMasses[7, 0], peptide.theoMasses[8, 0], peptide.theoMasses[9,0], peptide.theoMasses[10,0], peptide.theoMasses[11,0],
                            peptide.adjustedTheoMasses[0, 0], peptide.adjustedTheoMasses[1, 0], peptide.adjustedTheoMasses[2, 0], peptide.adjustedTheoMasses[3, 0], peptide.adjustedTheoMasses[4, 0], peptide.adjustedTheoMasses[5, 0], peptide.adjustedTheoMasses[6, 0], peptide.adjustedTheoMasses[7, 0], peptide.adjustedTheoMasses[8, 0], peptide.adjustedTheoMasses[9,0], peptide.adjustedTheoMasses[10,0], peptide.adjustedTheoMasses[11,0], peptide.GetTheoreticalResolvability(0), peptide.GetTheoreticalResolvability(1), peptide.GetTheoreticalResolvability(2),
                            peptide.countCompletePairs, peptide.countCompleteIsotopes[0], peptide.countCompleteIsotopes[1], peptide.countCompleteIsotopes[2], peptide.completeTotalIntensity[0, NUMISOTOPES], peptide.completeTotalIntensity[1, NUMISOTOPES], peptide.completeTotalIntensity[2, NUMISOTOPES], peptide.completeTotalIntensity[3, NUMISOTOPES], peptide.completeTotalIntensity[4, NUMISOTOPES], peptide.completeTotalIntensity[5, NUMISOTOPES], peptide.completeTotalIntensity[6, NUMISOTOPES], peptide.completeTotalIntensity[7, NUMISOTOPES], peptide.completeTotalIntensity[8, NUMISOTOPES], peptide.completeTotalIntensity[9,NUMISOTOPES], peptide.completeTotalIntensity[10,NUMISOTOPES], peptide.completeTotalIntensity[11, NUMISOTOPES],
                            peptide.heavyToLightRatioSum[0, 0], peptide.heavyToLightRatioSum[1, 0], peptide.heavyToLightRatioSum[2, 0], peptide.heavyToLightRatioSum[3, 0], peptide.heavyToLightRatioSum[4, 0], peptide.heavyToLightRatioSum[5, 0], peptide.heavyToLightRatioSum[6,0], peptide.heavyToLightRatioSum[7,0], peptide.heavyToLightRatioSum[8,0],peptide.finalQuantified[0], peptide.finalQuantified[1], peptide.finalQuantified[2], peptide.quantifiedNoiseIncluded[0], peptide.quantifiedNoiseIncluded[1], peptide.quantifiedNoiseIncluded[2]);
                    }
                }
                else
                {
                    header1 = ("Scan number(s), Raw File(s), Charge State(s), Best Scan Number, Peptide, # Labels, Theo Mass 1, Theo Mass 2, Theo Mass 3, Theo Mass 4, Theo Mass 5, Theo Mass 6, Theo Mass 7, Theo Mass 8, Theo Mass 9, Theo Mass 10, Theo Mass 11, Theo Mass 12, Adjusted Mass 1, Adjusted Mass 2, Adjusted Mass 3, Adjusted Mass 4, Adjusted Mass 5, Adjusted Mass 6, Adjusted Mass 7, Adjusted Mass 8, Adjusted Mass 9, Adjusted Mass 10, Adjusted Mass 11, Adjusted Mass 12, Cluster 1 Resolvable?, Cluster 2 Resolvable?, Cluster 3 Resolvable?, Quantified MS1 Scans, Cluster 1 # Measurements, Cluster 2 # Measurements, Cluster 3 # Measurements, Intensity 1, Intensity 2, Intensity 3, Intensity 4, Intensity 5, Intensity 6, Intensity 7, Intensity 8, Intensity 9, Intensity 10, Intensity 11, Intensity 12, Ratio 2/1, Ratio 3/1, Ratio 4/1, Ratio 6/5, Ratio 7/5, Ratio 8/5, Ratio 10/9, Ratio 11/9, Ratio 12/9, Cluster 1 Ratio Count, Cluster 2 Ratio Count, Cluster 3 Ratio Count, Cluster 1 Missing Channels?, Cluster 2 Missing Channels?, Cluster 3 Missing Channels?");
                    writer1.WriteLine(header1);
                    foreach (PeptideID peptide in allPeptides.Values)
                    {
                        writer1.WriteLine("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23},{24},{25},{26},{27},{28},{29},{30},{31},{32},{33},{34},{35},{36},{37},{38},{39},{40},{41},{42},{43},{44},{45},{46},{47},{48},{49},{50},{51},{52},{53},{54},{55},{56},{57},{58},{59},{60},{61},{62},{63}",
                            peptide.scanNumbers, peptide.rawFiles, peptide.chargeStates, peptide.bestPSM.ScanNumber, peptide.sequence, peptide.numLabels,
                            peptide.theoMasses[0, 0], peptide.theoMasses[1, 0], peptide.theoMasses[2, 0], peptide.theoMasses[3, 0], peptide.theoMasses[4, 0], peptide.theoMasses[5, 0], peptide.theoMasses[6, 0], peptide.theoMasses[7, 0], peptide.theoMasses[8, 0], peptide.theoMasses[9, 0], peptide.theoMasses[10, 0], peptide.theoMasses[11, 0],
                            peptide.adjustedTheoMasses[0, 0], peptide.adjustedTheoMasses[1, 0], peptide.adjustedTheoMasses[2, 0], peptide.adjustedTheoMasses[3, 0], peptide.adjustedTheoMasses[4, 0], peptide.adjustedTheoMasses[5, 0], peptide.adjustedTheoMasses[6, 0], peptide.adjustedTheoMasses[7, 0], peptide.adjustedTheoMasses[8, 0], peptide.adjustedTheoMasses[9, 0], peptide.adjustedTheoMasses[10, 0], peptide.adjustedTheoMasses[11, 0], peptide.GetTheoreticalResolvability(0), peptide.GetTheoreticalResolvability(1), peptide.GetTheoreticalResolvability(2),
                            peptide.countAllPairs, peptide.countAllIsotopes[0], peptide.countAllIsotopes[1], peptide.countAllIsotopes[2], peptide.totalIntensity[0, NUMISOTOPES], peptide.totalIntensity[1, NUMISOTOPES], peptide.totalIntensity[2, NUMISOTOPES], peptide.totalIntensity[3, NUMISOTOPES], peptide.totalIntensity[4, NUMISOTOPES], peptide.totalIntensity[5, NUMISOTOPES], peptide.totalIntensity[6, NUMISOTOPES], peptide.totalIntensity[7, NUMISOTOPES], peptide.totalIntensity[8, NUMISOTOPES], peptide.totalIntensity[9, NUMISOTOPES], peptide.totalIntensity[10, NUMISOTOPES], peptide.totalIntensity[11, NUMISOTOPES],
                            peptide.heavyToLightRatioSum[0, 0], peptide.heavyToLightRatioSum[1, 0], peptide.heavyToLightRatioSum[2, 0], peptide.heavyToLightRatioSum[3, 0], peptide.heavyToLightRatioSum[4, 0], peptide.heavyToLightRatioSum[5, 0], peptide.heavyToLightRatioSum[6, 0], peptide.heavyToLightRatioSum[7, 0], peptide.heavyToLightRatioSum[8, 0], peptide.finalQuantified[0], peptide.finalQuantified[1], peptide.finalQuantified[2], peptide.quantifiedNoiseIncluded[0], peptide.quantifiedNoiseIncluded[1], peptide.quantifiedNoiseIncluded[2]);
                    }
                }
                writer1.Close();
            }
            else if (NUMCHANNELS == 18)
            {
                string header1;
                if (!NOISEBANDCAP)
                {
                    header1 = ("Scan number(s), Raw File(s), Charge State(s), Best Scan Number, Peptide, # Labels, Theo Mass 1, Theo Mass 2, Theo Mass 3, Theo Mass 4, Theo Mass 5, Theo Mass 6, Theo Mass 7, Theo Mass 8, Theo Mass 9, Theo Mass 10, Theo Mass 11, Theo Mass 12, Theo Mass 13, Theo Mass 14, Theo Mass 15, Theo Mass 16, Theo Mass 17, Theo  Mass 18, Adjusted Mass 1, Adjusted Mass 2, Adjusted Mass 3, Adjusted Mass 4, Adjusted Mass 5, Adjusted Mass 6, Adjusted Mass 7, Adjusted Mass 8, Adjusted Mass 9, Adjusted Mass 10, Adjusted Mass 11, Adjusted Mass 12, Adjusted Mass 13, Adjusted Mass 14, Adjusted Mass 15, Adjusted Mass 16, Adjusted Mass 17, Adjusted Mass 18, Cluster 1 Resolvable?, Cluster 2 Resolvable?, Cluster 3 Resolvable?, Quantified MS1 Scans, Cluster 1 # Measurements, Cluster 2 # Measurements, Cluster 3 # Measurements, Intensity 1, Intensity 2, Intensity 3, Intensity 4, Intensity 5, Intensity 6, Intensity 7, Intensity 8, Intensity 9, Intensity 10, Intensity 11, Intensity 12, Intensity 13, Intensity 14, Intensity 15, Intensity 16, Intensity 17, Intensity 18, Ratio 2/1, Ratio 3/1, Ratio 4/1, Ratio 5/1, Ratio 6/1, Ratio 8/7, Ratio 9/7, Ratio 10/7, Ratio 11/7, Ratio 12/7, Ratio 14/13, Ratio 15/13, Ratio 16/13, Ratio 17/13, Ratio 18/13, Cluster 1 Ratio Count, Cluster 2 Ratio Count, Cluster 3 Ratio Count, Cluster 1 Missing Channels?, Cluster 2 Missing Channels?, Cluster 3 Missing Channels?");
                    writer1.WriteLine(header1);
                    foreach (PeptideID peptide in allPeptides.Values)
                    {
                        writer1.WriteLine("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23},{24},{25},{26},{27},{28},{29},{30},{31},{32},{33},{34},{35},{36},{37},{38},{39},{40},{41},{42},{43},{44},{45},{46},{47},{48},{49},{50},{51},{52},{53},{54},{55},{56},{57},{58},{59},{60},{61},{62},{63},{64},{65},{66},{67},{68},{69},{70},{71},{72},{73},{74},{75},{76},{77},{78},{79},{80},{81},{82},{83},{84},{85},{86},{87}",
                            peptide.scanNumbers, peptide.rawFiles, peptide.chargeStates, peptide.bestPSM.ScanNumber, peptide.sequence, peptide.numLabels,
                            peptide.theoMasses[0, 0], peptide.theoMasses[1, 0], peptide.theoMasses[2, 0], peptide.theoMasses[3, 0], peptide.theoMasses[4, 0], peptide.theoMasses[5, 0], peptide.theoMasses[6, 0], peptide.theoMasses[7, 0], peptide.theoMasses[8, 0], peptide.theoMasses[9, 0], peptide.theoMasses[10, 0], peptide.theoMasses[11, 0], peptide.theoMasses[12,0], peptide.theoMasses[13,0], peptide.theoMasses[14,0], peptide.theoMasses[15,0], peptide.theoMasses[16,0], peptide.theoMasses[17,0],
                            peptide.adjustedTheoMasses[0, 0], peptide.adjustedTheoMasses[1, 0], peptide.adjustedTheoMasses[2, 0], peptide.adjustedTheoMasses[3, 0], peptide.adjustedTheoMasses[4, 0], peptide.adjustedTheoMasses[5, 0], peptide.adjustedTheoMasses[6, 0], peptide.adjustedTheoMasses[7, 0], peptide.adjustedTheoMasses[8, 0], peptide.adjustedTheoMasses[9, 0], peptide.adjustedTheoMasses[10, 0], peptide.adjustedTheoMasses[11, 0], peptide.adjustedTheoMasses[12,0], peptide.adjustedTheoMasses[13,0], peptide.adjustedTheoMasses[14,0], peptide.adjustedTheoMasses[15,0], peptide.adjustedTheoMasses[16,0], peptide.adjustedTheoMasses[17,0], peptide.GetTheoreticalResolvability(0), peptide.GetTheoreticalResolvability(1), peptide.GetTheoreticalResolvability(2),
                            peptide.countCompletePairs, peptide.countCompleteIsotopes[0], peptide.countCompleteIsotopes[1], peptide.countCompleteIsotopes[2], peptide.completeTotalIntensity[0, NUMISOTOPES], peptide.completeTotalIntensity[1, NUMISOTOPES], peptide.completeTotalIntensity[2, NUMISOTOPES], peptide.completeTotalIntensity[3, NUMISOTOPES], peptide.completeTotalIntensity[4, NUMISOTOPES], peptide.completeTotalIntensity[5, NUMISOTOPES], peptide.completeTotalIntensity[6, NUMISOTOPES], peptide.completeTotalIntensity[7, NUMISOTOPES], peptide.completeTotalIntensity[8, NUMISOTOPES], peptide.completeTotalIntensity[9, NUMISOTOPES], peptide.completeTotalIntensity[10, NUMISOTOPES], peptide.completeTotalIntensity[11, NUMISOTOPES], peptide.completeTotalIntensity[12, NUMISOTOPES], peptide.completeTotalIntensity[13, NUMISOTOPES], peptide.completeTotalIntensity[14,0], peptide.completeTotalIntensity[15,0], peptide.completeTotalIntensity[16,0], peptide.completeTotalIntensity[17,0],
                            peptide.heavyToLightRatioSum[0, 0], peptide.heavyToLightRatioSum[1, 0], peptide.heavyToLightRatioSum[2, 0], peptide.heavyToLightRatioSum[3, 0], peptide.heavyToLightRatioSum[4, 0], peptide.heavyToLightRatioSum[5, 0], peptide.heavyToLightRatioSum[6, 0], peptide.heavyToLightRatioSum[7, 0], peptide.heavyToLightRatioSum[8, 0], peptide.heavyToLightRatioSum[9,0], peptide.heavyToLightRatioSum[10,0], peptide.heavyToLightRatioSum[11,0], peptide.heavyToLightRatioSum[12,0], peptide.heavyToLightRatioSum[13,0], peptide.heavyToLightRatioSum[14,0], peptide.finalQuantified[0], peptide.finalQuantified[1], peptide.finalQuantified[2], peptide.quantifiedNoiseIncluded[0], peptide.quantifiedNoiseIncluded[1], peptide.quantifiedNoiseIncluded[2]);
                    }
                }
                else
                {
                    header1 = ("Scan number(s), Raw File(s), Charge State(s), Best Scan Number, Peptide, # Labels, Theo Mass 1, Theo Mass 2, Theo Mass 3, Theo Mass 4, Theo Mass 5, Theo Mass 6, Theo Mass 7, Theo Mass 8, Theo Mass 9, Theo Mass 10, Theo Mass 11, Theo Mass 12, Theo Mass 13, Theo Mass 14, Theo Mass 15, Theo Mass 16, Theo Mass 17, Theo  Mass 18, Adjusted Mass 1, Adjusted Mass 2, Adjusted Mass 3, Adjusted Mass 4, Adjusted Mass 5, Adjusted Mass 6, Adjusted Mass 7, Adjusted Mass 8, Adjusted Mass 9, Adjusted Mass 10, Adjusted Mass 11, Adjusted Mass 12, Adjusted Mass 13, Adjusted Mass 14, Adjusted Mass 15, Adjusted Mass 16, Adjusted Mass 17, Adjusted Mass 18, Cluster 1 Resolvable?, Cluster 2 Resolvable?, Cluster 3 Resolvable?, Quantified MS1 Scans, Cluster 1 # Measurements, Cluster 2 # Measurements, Cluster 3 # Measurements, Intensity 1, Intensity 2, Intensity 3, Intensity 4, Intensity 5, Intensity 6, Intensity 7, Intensity 8, Intensity 9, Intensity 10, Intensity 11, Intensity 12, Intensity 13, Intensity 14, Intensity 15, Intensity 16, Intensity 17, Intensity 18, Ratio 2/1, Ratio 3/1, Ratio 4/1, Ratio 5/1, Ratio 6/1, Ratio 8/7, Ratio 9/7, Ratio 10/7, Ratio 11/7, Ratio 12/7, Ratio 14/13, Ratio 15/13, Ratio 16/13, Ratio 17/13, Ratio 18/13, Cluster 1 Ratio Count, Cluster 2 Ratio Count, Cluster 3 Ratio Count, Cluster 1 Missing Channels?, Cluster 2 Missing Channels?, Cluster 3 Missing Channels?");
                    writer1.WriteLine(header1);
                    foreach (PeptideID peptide in allPeptides.Values)
                    {
                        writer1.WriteLine("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23},{24},{25},{26},{27},{28},{29},{30},{31},{32},{33},{34},{35},{36},{37},{38},{39},{40},{41},{42},{43},{44},{45},{46},{47},{48},{49},{50},{51},{52},{53},{54},{55},{56},{57},{58},{59},{60},{61},{62},{63},{64},{65},{66},{67},{68},{69},{70},{71},{72},{73},{74},{75},{76},{77},{78},{79},{80},{81},{82},{83},{84},{85},{86},{87}",
                            peptide.scanNumbers, peptide.rawFiles, peptide.chargeStates, peptide.bestPSM.ScanNumber, peptide.sequence, peptide.numLabels,
                            peptide.theoMasses[0, 0], peptide.theoMasses[1, 0], peptide.theoMasses[2, 0], peptide.theoMasses[3, 0], peptide.theoMasses[4, 0], peptide.theoMasses[5, 0], peptide.theoMasses[6, 0], peptide.theoMasses[7, 0], peptide.theoMasses[8, 0], peptide.theoMasses[9, 0], peptide.theoMasses[10, 0], peptide.theoMasses[11, 0], peptide.theoMasses[12, 0], peptide.theoMasses[13, 0], peptide.theoMasses[14, 0], peptide.theoMasses[15, 0], peptide.theoMasses[16, 0], peptide.theoMasses[17, 0],
                            peptide.adjustedTheoMasses[0, 0], peptide.adjustedTheoMasses[1, 0], peptide.adjustedTheoMasses[2, 0], peptide.adjustedTheoMasses[3, 0], peptide.adjustedTheoMasses[4, 0], peptide.adjustedTheoMasses[5, 0], peptide.adjustedTheoMasses[6, 0], peptide.adjustedTheoMasses[7, 0], peptide.adjustedTheoMasses[8, 0], peptide.adjustedTheoMasses[9, 0], peptide.adjustedTheoMasses[10, 0], peptide.adjustedTheoMasses[11, 0], peptide.adjustedTheoMasses[12, 0], peptide.adjustedTheoMasses[13, 0], peptide.adjustedTheoMasses[14, 0], peptide.adjustedTheoMasses[15, 0], peptide.adjustedTheoMasses[16, 0], peptide.adjustedTheoMasses[17, 0], peptide.GetTheoreticalResolvability(0), peptide.GetTheoreticalResolvability(1), peptide.GetTheoreticalResolvability(2),
                            peptide.countAllPairs, peptide.countAllIsotopes[0], peptide.countAllIsotopes[1], peptide.countAllIsotopes[2], peptide.totalIntensity[0, NUMISOTOPES], peptide.totalIntensity[1, NUMISOTOPES], peptide.totalIntensity[2, NUMISOTOPES], peptide.totalIntensity[3, NUMISOTOPES], peptide.totalIntensity[4, NUMISOTOPES], peptide.totalIntensity[5, NUMISOTOPES], peptide.totalIntensity[6, NUMISOTOPES], peptide.totalIntensity[7, NUMISOTOPES], peptide.totalIntensity[8, NUMISOTOPES], peptide.totalIntensity[9, NUMISOTOPES], peptide.totalIntensity[10, NUMISOTOPES], peptide.totalIntensity[11, NUMISOTOPES], peptide.totalIntensity[12, NUMISOTOPES], peptide.totalIntensity[13, NUMISOTOPES], peptide.totalIntensity[14, 0], peptide.totalIntensity[15, 0], peptide.totalIntensity[16, 0], peptide.totalIntensity[17, 0],
                            peptide.heavyToLightRatioSum[0, 0], peptide.heavyToLightRatioSum[1, 0], peptide.heavyToLightRatioSum[2, 0], peptide.heavyToLightRatioSum[3, 0], peptide.heavyToLightRatioSum[4, 0], peptide.heavyToLightRatioSum[5, 0], peptide.heavyToLightRatioSum[6, 0], peptide.heavyToLightRatioSum[7, 0], peptide.heavyToLightRatioSum[8, 0], peptide.heavyToLightRatioSum[9, 0], peptide.heavyToLightRatioSum[10, 0], peptide.heavyToLightRatioSum[11, 0], peptide.heavyToLightRatioSum[12, 0], peptide.heavyToLightRatioSum[13, 0], peptide.heavyToLightRatioSum[14, 0], peptide.finalQuantified[0], peptide.finalQuantified[1], peptide.finalQuantified[2], peptide.quantifiedNoiseIncluded[0], peptide.quantifiedNoiseIncluded[1], peptide.quantifiedNoiseIncluded[2]);
                    }
                }
                writer1.Close();
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
            //else if (NEUCODE_SIXPLEX_MTRAQ)
            //{
            //    string header2 = ("Scan number(s), Raw File(s), Charge State(s), Best Scan Number, Peptide, # Labels, Theo Mass 1, Theo Mass 2, Theo Mass 3, Theo Mass 4, Theo Mass 5, Theo Mass 6, Adjusted Mass 1, Adjusted Mass 2, Adjusted Mass 3, Adjusted Mass 4, Adjusted Mass 5, Adjusted Mass 6, Quantified MS1 Scans, # Total Measurements (Rep 1), # Total Measurements (Rep 2), # Total Measurements (Rep 3), Total Intensity 1, Total Intensity 2, Total Intensity 3, Total Intensity 4, Total Intensity 5, Total Intensity 6, Quantified no NBC MS1 Scans, # No NBC Measurements (Rep 1), # No NBC Measurements (Rep 2), # No NBC Measurements (Rep 3), No NBC Intensity 1, No NBC Intensity 2, No NBC Intensity 3, No NBC Intensity 4, No NBC Intensity 5, No NBC Intensity 6, Coalescence Detected?, Ratio (Rep 1), Ratio Count (Rep 1), Missing Channels Quantified (Rep 1), Ratio (Rep 2), Ratio Count (Rep 2), Missing Channels Quantified (Rep 2), Ratio (Rep 3), Ratio Count (Rep 3), Missing Channels Quantified (Rep 3)");
            //    writer1.WriteLine(header2);

            //    foreach (PeptideID peptide in allPeptides.Values)
            //    {
            //        writer1.WriteLine("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23},{24},{25},{26},{27},{28},{29},{30},{31},{32},{33},{34},{35},{36},{37},{38},{39},{40},{41},{42},{43},{44},{45},{46},{47}",
            //            peptide.scanNumbers, peptide.rawFiles, peptide.chargeStates, peptide.bestPSM.ScanNumber, peptide.sequence, peptide.numLabels,
            //            peptide.theoMasses[0, 0], peptide.theoMasses[1, 0], peptide.theoMasses[2, 0], peptide.theoMasses[3, 0], peptide.theoMasses[4, 0], peptide.theoMasses[5, 0],
            //            peptide.adjustedTheoMasses[0, 0], peptide.adjustedTheoMasses[1, 0], peptide.adjustedTheoMasses[2, 0], peptide.adjustedTheoMasses[3, 0], peptide.adjustedTheoMasses[4, 0], peptide.adjustedTheoMasses[5, 0],
            //            peptide.countAllPairs, peptide.countAllIsotopes[0], peptide.countAllIsotopes[1], peptide.countAllIsotopes[2],
            //            peptide.totalIntensity[0, NUMISOTOPES], peptide.totalIntensity[1, NUMISOTOPES], peptide.totalIntensity[2, NUMISOTOPES], peptide.totalIntensity[3, NUMISOTOPES], peptide.totalIntensity[4, NUMISOTOPES], peptide.totalIntensity[5, NUMISOTOPES],
            //            peptide.countCompletePairs, peptide.countCompleteIsotopes[0], peptide.countCompleteIsotopes[1], peptide.countCompleteIsotopes[2],
            //            peptide.completeTotalIntensity[0, NUMISOTOPES], peptide.completeTotalIntensity[1, NUMISOTOPES], peptide.completeTotalIntensity[2, NUMISOTOPES], peptide.completeTotalIntensity[3, NUMISOTOPES], peptide.completeTotalIntensity[4, NUMISOTOPES], peptide.completeTotalIntensity[5, NUMISOTOPES],
            //            peptide.coalescenceDetected,
            //            peptide.heavyToLightRatioSum[0, 0], peptide.finalQuantified[0], peptide.quantifiedNoiseIncluded[0],
            //            peptide.heavyToLightRatioSum[1, 0], peptide.finalQuantified[1], peptide.quantifiedNoiseIncluded[1],
            //            peptide.heavyToLightRatioSum[2, 0], peptide.finalQuantified[2], peptide.quantifiedNoiseIncluded[2]);
            //    }
            //}
            //else
            //{
            //    string header3 = ("Scan number(s), Raw File(s), Charge State(s), Best Scan Number, Peptide, # K, # R, Light R Light Theo Mass, Light R Heavy Theo Mass, Light R Adjusted Light Theo Mass, Light R Adjusted Heavy Theo Mass, Medium R Light Theo Mass, Medium R Heavy Theo Mass, Medium R Adjusted Light Theo Mass, Medium R Adjusted Heavy Theo Mass, Heavy R Light Theo Mass, Heavy R Heavy Theo Mass, Heavy R Adjusted Light Theo Mass, Heavy R Adjusted Heavy Theo Mass, Quantified MS1 Scans, # Total Measurements, Total Light R Light Intensity, Total Light R Heavy Intensity, Total Medium R Light Intensity, Total Medium R Heavy Intensity, Total Heavy R Light Intensity, Total Heavy R Heavy Intensity, Quantified no NBC MS1 Scans, # No NBC Measurements (Rep 1), # No NBC Measurements (Rep 2), # No NBC Measurements (Rep 3), Light R Light no NBC Intensity, Light R Heavy no NBC Intensity, Medium R Light no NBC Intensity, Medium R Heavy no NBC Intensity, Heavy R Light no NBC Intensity, Coalescence Detected?, Ratio (Rep 1), Ratio Count (Rep 1), Missing Channels Quantified (Rep 1), Ratio (Rep 2), Ratio Count (Rep 2), Missing Channels Quantified (Rep 2), Ratio (Rep 3), Ratio Count (Rep 3), Missing Channels Quantified (Rep 3)");
            //    writer1.WriteLine(header3);

            //    foreach (PeptideID peptide in allPeptides.Values)
            //    {
            //        if (peptide.countResidues('R', peptide.sequence) == 0)
            //        {
            //            writer1.WriteLine("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23},{24},{25},{26},{27},{28},{29},{30},{31},{32},{33},{34},{35},{36},{37},{38},{39},{40}",
            //            peptide.scanNumbers, peptide.rawFiles, peptide.chargeStates, peptide.bestPSM.ScanNumber, peptide.sequence, peptide.countResidues('K', peptide.sequence), 0,
            //            peptide.theoMasses[0, 0], peptide.theoMasses[1, 0], 0, 0, 0, 0,
            //            peptide.adjustedTheoMasses[0, 0], peptide.adjustedTheoMasses[1, 0], 0, 0, 0, 0,
            //            peptide.countAllPairs, peptide.countAllIsotopes, peptide.totalIntensity[0, NUMISOTOPES], peptide.totalIntensity[1, NUMISOTOPES], 0, 0, 0, 0,
            //            peptide.countCompletePairs, peptide.countCompleteIsotopes, peptide.completeTotalIntensity[0, NUMISOTOPES], peptide.completeTotalIntensity[1, NUMISOTOPES], 0, 0, 0, 0,
            //            peptide.coalescenceDetected, peptide.heavyToLightRatioSum[0, 0], 0, 0, peptide.finalQuantified, peptide.quantifiedNoiseIncluded);
            //        }
            //        else
            //        {
            //            writer1.WriteLine("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23},{24},{25},{26},{27},{28},{29},{30},{31},{32},{33},{34},{35},{36},{37},{38},{39},{40}",
            //            peptide.scanNumbers, peptide.rawFiles, peptide.chargeStates, peptide.bestPSM.ScanNumber, peptide.sequence, peptide.countResidues('K', peptide.sequence), peptide.countResidues('R', peptide.sequence),
            //            peptide.theoMasses[0, 0], peptide.theoMasses[1, 0], peptide.theoMasses[2, 0], peptide.theoMasses[3, 0], peptide.theoMasses[4, 0], peptide.theoMasses[5, 0],
            //            peptide.adjustedTheoMasses[0, 0], peptide.adjustedTheoMasses[1, 0], peptide.adjustedTheoMasses[2, 0], peptide.adjustedTheoMasses[3, 0], peptide.adjustedTheoMasses[4, 0], peptide.adjustedTheoMasses[5, 0],
            //            peptide.countAllPairs, peptide.countAllIsotopes, peptide.totalIntensity[0, NUMISOTOPES], peptide.totalIntensity[1, NUMISOTOPES], peptide.totalIntensity[2, NUMISOTOPES], peptide.totalIntensity[3, NUMISOTOPES], peptide.totalIntensity[4, NUMISOTOPES], peptide.totalIntensity[5, NUMISOTOPES],
            //            peptide.countCompletePairs, peptide.countCompleteIsotopes, peptide.completeTotalIntensity[0, NUMISOTOPES], peptide.completeTotalIntensity[1, NUMISOTOPES], peptide.completeTotalIntensity[2, NUMISOTOPES], peptide.completeTotalIntensity[3, NUMISOTOPES], peptide.completeTotalIntensity[4, NUMISOTOPES], peptide.completeTotalIntensity[5, NUMISOTOPES],
            //            peptide.coalescenceDetected, peptide.heavyToLightRatioSum[0, 0], peptide.heavyToLightRatioSum[1, 0], peptide.heavyToLightRatioSum[2, 0], peptide.finalQuantified, peptide.quantifiedNoiseIncluded);
            //        }
            //    }
            //}
            writer1.Close();
        }

        /* Reads in peptide information from a CSV file
         * Necessary information: spectrum number, charge, peptide sequence, E-value, raw file
         */
        private void readCsvInputFile(Dictionary<string, PeptideID> allPeptides, string file)
        {
            // Amino acids
            NamedChemicalFormula K600 = NamedChemicalFormula.AddModification("C-6 C{13}6", "Lys +6 13C6");
            NamedChemicalFormula K602 = NamedChemicalFormula.AddModification("C-6 C{13}6 N-2 N{15}2", "Lys +8 13C6 15N2");
            NamedChemicalFormula K422 = NamedChemicalFormula.AddModification("C-4 C{13}4 H-2 H{2}2 N-2 N{15}2", "Lys +8 13C4 2H2 15N2");
            NamedChemicalFormula K521 = NamedChemicalFormula.AddModification("C-5 C{13}5 H-2 H{2}2 N-1 N{15}1", "Lys +8 13C5 2H2 15N1");
            NamedChemicalFormula K062 = NamedChemicalFormula.AddModification("H-6 H{2}6 N-2 N{15}2", "Lys +8 2H6 15N2");
            NamedChemicalFormula K440 = NamedChemicalFormula.AddModification("C-4 C{13}4 H-4 H{2}4", "Lys +8 13C4 2H4");
            NamedChemicalFormula K080 = NamedChemicalFormula.AddModification("H-8 H{2}8", "Lys +8 2H8");
            NamedChemicalFormula R600 = NamedChemicalFormula.AddModification("C-6 C{13}6", "Arg +6 13C6");
            NamedChemicalFormula R604 = NamedChemicalFormula.AddModification("C-6 C{13}6 N-4 N{15}4", "Arg +10 13C6 15N4");
            NamedChemicalFormula K100 = NamedChemicalFormula.AddModification("C-1 C{13}", "Lys +1 13C");
            NamedChemicalFormula K001 = NamedChemicalFormula.AddModification("N-1 N{15}", "Lys +1 15N");
            NamedChemicalFormula L601 = NamedChemicalFormula.AddModification("N-1 N{15}1 C-6 C{13}6", "Leu +7 13C6 15N1");
            NamedChemicalFormula L070 = NamedChemicalFormula.AddModification("H-7 H{2}7", "Leu +7 2H7");
            NamedChemicalFormula L600 = NamedChemicalFormula.AddModification("C-6 C{13}6", "Leu +6 13C6");
            NamedChemicalFormula L0100 = NamedChemicalFormula.AddModification("H-10 H{2}10", "Leu +10 2H10");

            // mTRAQ labels
            NamedChemicalFormula lightmTRAQ = NamedChemicalFormula.AddModification("H{1}12 C{12}7 N{14}2 O{16}1", "mTRAQ L");
            NamedChemicalFormula lightmTRAQK602 = NamedChemicalFormula.AddModification("H{1}12 C{12}7 N{14}2 O{16}1 C-6 C{13}6 N-2 N{15}2", "mTRAQ L Lys +8 13C6 15N2");
            NamedChemicalFormula lightmTRAQK422 = NamedChemicalFormula.AddModification("H{1}12 C{12}7 N{14}2 O{16}1 C-4 C{13}4 H-4 H{2}4 N-2 N{15}2", "mTRAQ L Lys +8 13C4 2H2 15N2");
            NamedChemicalFormula lightmTRAQK521 = NamedChemicalFormula.AddModification("H{1}12 C{12}7 N{14}2 O{16}1 C-6 C{13}5 H-2 H{2}2 N-1 N{15}1", "mTRAQ L Lys +8 13C5 2H2 15N1");
            NamedChemicalFormula lightmTRAQK062 = NamedChemicalFormula.AddModification("H{1}12 C{12}7 N{14}2 O{16}1 H-6 H{2}6 N-2 N{15}2", "mTRAQ L Lys +8 2H6 15N2");
            NamedChemicalFormula lightmTRAQK440 = NamedChemicalFormula.AddModification("H{1}12 C{12}7 N{14}2 O{16}1 C-4 C{13}4 H-4 H{2}4", "mTRAQ L Lys +8 13C4 2H4");
            NamedChemicalFormula lightmTRAQK080 = NamedChemicalFormula.AddModification("H{1}12 C{12}7 N{14}2 O{16}1 H-8 H{2}8", "mTRAQ L Lys +8 2H8");
            NamedChemicalFormula mediummTRAQ = NamedChemicalFormula.AddModification("H{1}12 C{12}4 C{13}3 N{14}1 N{15}1 O{16}1", "mTRAQ M");
            NamedChemicalFormula mediummTRAQ602 = NamedChemicalFormula.AddModification("H{1}12 C{12}4 C{13}3 N{14}1 N{15}1 O{16}1 C-6 C{13}6 N-2 N{15}2", "mTRAQ M Lys +8 13C6 15N2");
            NamedChemicalFormula mediummTRAQK422 = NamedChemicalFormula.AddModification("H{1}12 C{12}4 C{13}3 N{14}1 N{15}1 O{16}1 C-4 C{13}4 H-4 H{2}4 N-2 N{15}2", "mTRAQ M Lys +8 13C4 2H2 15N2");
            NamedChemicalFormula mediummTRAQK521 = NamedChemicalFormula.AddModification("H{1}12 C{12}4 C{13}3 N{14}1 N{15}1 O{16}1 C-6 C{13}5 H-2 H{2}2 N-1 N{15}1", "mTRAQ M Lys +8 13C5 2H2 15N1");
            NamedChemicalFormula mediummTRAQK062 = NamedChemicalFormula.AddModification("H{1}12 C{12}4 C{13}3 N{14}1 N{15}1 O{16}1 H-6 H{2}6 N-2 N{15}2", "mTRAQ M Lys +8 2H6 15N2");
            NamedChemicalFormula mediummTRAQK440 = NamedChemicalFormula.AddModification("H{1}12 C{12}4 C{13}3 N{14}1 N{15}1 O{16}1 C-4 C{13}4 H-4 H{2}4", "mTRAQ M Lys +8 13C4 2H4");
            NamedChemicalFormula mediummTRAQK080 = NamedChemicalFormula.AddModification("H{1}12 C{12}4 C{13}3 N{14}1 N{15}1 O{16}1 H-8 H{2}8", "mTRAQ M Lys +8 2H8");
            NamedChemicalFormula heavymTRAQ = NamedChemicalFormula.AddModification("H{1}12 C{12}1 C{13}6 N{15}2 O{16}1", "mTRAQ H");
            NamedChemicalFormula heavymTRAQK602 = NamedChemicalFormula.AddModification("H{1}12 C{12}1 C{13}6 N{15}2 O{16}1 C-6 C{13}6 N-2 N{15}2", "mTRAQ H Lys +8 13C6 15N2");
            NamedChemicalFormula heavymTRAQK422 = NamedChemicalFormula.AddModification("H{1}12 C{12}1 C{13}6 N{15}2 O{16}1 C-4 C{13}4 H-4 H{2}4 N-2 N{15}2", "mTRAQ H Lys +8 13C4 2H2 15N2");
            NamedChemicalFormula heavymTRAQK521 = NamedChemicalFormula.AddModification("H{1}12 C{12}1 C{13}6 N{15}2 O{16}1 C-6 C{13}5 H-2 H{2}2 N-1 N{15}1", "mTRAQ H Lys +8 13C5 2H2 15N1");
            NamedChemicalFormula heavymTRAQK062 = NamedChemicalFormula.AddModification("H{1}12 C{12}1 C{13}6 N{15}2 O{16}1 H-6 H{2}6 N-2 N{15}2", "mTRAQ H Lys +8 2H6 15N2");
            NamedChemicalFormula heavymTRAQK440 = NamedChemicalFormula.AddModification("H{1}12 C{12}1 C{13}6 N{15}2 O{16}1 C-4 C{13}4 H-4 H{2}4", "mTRAQ H Lys +8 13C4 2H4");
            NamedChemicalFormula heavymTRAQK080 = NamedChemicalFormula.AddModification("H{1}12 C{12}1 C{13}6 N{15}2 O{16}1 H-8 H{2}8", "mTRAQ H Lys +8 2H8");
            
            // Chemical labels
            NamedChemicalFormula lightASH1 = NamedChemicalFormula.AddModification("C{12}18 H{1}31 N{14}1 O{16}5 N{15}6", "4plex L1");
            NamedChemicalFormula lightASH2 = NamedChemicalFormula.AddModification("C{12}16 H{1}31 N{14}3 O{16}5 N{15}4 C{13}2", "4plex L2");
            NamedChemicalFormula lightASH3 = NamedChemicalFormula.AddModification("C{12}14 H{1}31 N{14}5 O{16}5 N{15}2 C{13}4", "4plex L3");
            NamedChemicalFormula lightASH4 = NamedChemicalFormula.AddModification("C{12}12 H{1}31 N{14}7 O{16}5 C{13}6", "4plex L4");
            NamedChemicalFormula mediumASH1 = NamedChemicalFormula.AddModification("C{12}14 H{1}31 N{14}1 O{16}5 N{15}6 C{13}4", "4plex M1");
            NamedChemicalFormula mediumASH2 = NamedChemicalFormula.AddModification("C{12}12 H{1}31 N{14}1 O{16}5 N{15}4 C{13}6", "4plex M2");
            NamedChemicalFormula mediumASH3 = NamedChemicalFormula.AddModification("C{12}10 H{1}31 N{14}1 O{16}5 N{15}2 C{13}8", "4plex M3");
            NamedChemicalFormula mediumASH4 = NamedChemicalFormula.AddModification("C{12}8 H{1}31 N{14}1 O{16}5 C{13}10", "4plex M4");
            NamedChemicalFormula heavyASH1 = NamedChemicalFormula.AddModification("C{12}10 H{1}31 N{14}1 O{16}5 N{15}6 C{13}8", "4plex H1");
            NamedChemicalFormula heavyASH2 = NamedChemicalFormula.AddModification("C{12}8 H{1}31 N{14}3 O{16}5 N{15}4 C{13}10", "4plex H2");
            NamedChemicalFormula heavyASH3 = NamedChemicalFormula.AddModification("C{12}6 H{1}31 N{14}5 O{16}5 N{15}2 C{13}12", "4plex H3");
            NamedChemicalFormula heavyASH4 = NamedChemicalFormula.AddModification("C{12}4 H{1}31 N{14}7 O{16}5 C{13}14", "4plex H4");
            NamedChemicalFormula lightCarbamyl = NamedChemicalFormula.AddModification("C{12}1 N{15}1 H{1}2 O{16}1", "Carbamyl L");
            NamedChemicalFormula heavyCarbamyl = NamedChemicalFormula.AddModification("C{13}1 N{14}1 H{1}2 O{16}1", "Carbamyl H");
            NamedChemicalFormula lightICAT = NamedChemicalFormula.AddModification("C{12}10H{1}15N{14}3O{16}2", "iCAT L");
            NamedChemicalFormula heavyICAT = NamedChemicalFormula.AddModification("C{13}9C{12}1H{1}15N{14}3O{16}2", "iCAT H");

            HashSet<string> rawFiles = new HashSet<string>();

            //Cycle through .csv file to make a list of identified peptides and properties
            CsvReader targetFileReader = new CsvReader(new StreamReader(file), true);
            CsvReader injectionTimesReader;

            using (targetFileReader)
            {
                PeptideID newPeptide;
                MSDataFile rawFile = null;
                WriteMessage("uploading peptide IDs");
                while (targetFileReader.ReadNextRecord())
                {
                    string basePathName = targetFileReader["Filename/id"].Substring(0, targetFileReader["Filename/id"].IndexOf("."));

                    if (!RAWFILES.TryGetValue(basePathName, out rawFile))
                    {
                        rawFile = new ThermoRawFile(Path.Combine(rawFileBox.Text, basePathName + ".raw"));
                        RAWFILES.Add(basePathName, rawFile);

                        if (SEGMENTEDINJECTIONTIMES)
                        {
                            string injectionTimesFile = (Path.Combine(rawFileBox.Text, basePathName + "_times.csv"));
                            injectionTimesReader = new CsvReader(new StreamReader(injectionTimesFile), true);
                            INJECTIONTIMES.Add(basePathName, new Dictionary<int, Dictionary<Range<double>, double>>());

                            int scanNumber;
                            double firstMass;
                            double lastMass;
                            double injectionTime;

                            while (injectionTimesReader.ReadNextRecord())
                            {
                                Dictionary<int, Dictionary<Range<double>, double>> scanNumbers = null;
                                
                                scanNumber = int.Parse(injectionTimesReader["Scan Number"]);
                                firstMass = double.Parse(injectionTimesReader["First Mass"]);
                                lastMass = double.Parse(injectionTimesReader["Last Mass"]);
                                injectionTime = double.Parse(injectionTimesReader["Time"]);

                                INJECTIONTIMES.TryGetValue(basePathName, out scanNumbers);
                                Dictionary<Range<double>, double> segmentInjectionTimes = null;

                                if (scanNumbers.TryGetValue(scanNumber, out segmentInjectionTimes))
                                {
                                    segmentInjectionTimes.Add(new Range<double>(firstMass, lastMass), injectionTime);
                                }
                                else
                                {
                                    segmentInjectionTimes = new Dictionary<Range<double>, double>();
                                    segmentInjectionTimes.Add(new Range<double>(firstMass, lastMass), injectionTime);
                                    scanNumbers.Add(scanNumber, segmentInjectionTimes);
                                }
                            }
                        }
                    }

                    newPeptide = new PeptideID(int.Parse(targetFileReader["Spectrum number"]),
                                int.Parse(targetFileReader["Charge"]), double.Parse(targetFileReader["E-value"]),
                                targetFileReader["Peptide"], rawFile, targetFileReader["Mods"]);
                    checkAdd(allPeptides, newPeptide, (MSDataFile)rawFile, double.Parse(targetFileReader["E-value"]));                 
                }
                WriteMessage("done uploading peptide IDs: " + allPeptides.Count + " total unique sequences");
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
                    if (peptide.PSM.EValue < peptideID.bestPSMs[rawFile].EValue)
                    {
                        peptideID.bestPSMs[rawFile] = peptide.PSM;
                    }
                }
                else 
                {
                    // Raw file not in PSM dictionary
                    peptideID.PSMs.Add(rawFile, new List<PeptideSpectralMatch>());
                    List<PeptideSpectralMatch> psms = null;
                    peptideID.PSMs.TryGetValue(rawFile, out psms);
                    psms.Add(peptide.PSM);
                    peptideID.bestPSMs.Add(rawFile, peptide.PSM);
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
            WriteMessage("writing output");*/
    
            //Retrieve each peptide's maximum intensity & missing channel frequency
            foreach (PeptideID uniquePeptide in allPeptides.Values)
            {
                double maximumIntensity = uniquePeptide.log10MaxIntensity;
                double missingChannelFrequency = uniquePeptide.missingChannelFrequency;
                
                if (uniquePeptide.log10MaxIntensity > 0 && missingChannelFrequency >= 0)
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
                listBox1.Items.AddRange(browseTargetInput.FileNames);
                //csvInputBox.Text = browseTargetInput.FileName;
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
            start.Enabled = false;
            foreach (string file in listBox1.Items)
            {
                Run(file);
            }
            start.Enabled = true;            
        }

        private void Form1_Load(object sender, EventArgs e)
        {

        }

        private void Form1_DragEnter(object sender, DragEventArgs e)
        {
            e.Effect = DragDropEffects.All;
        }

        private void Form1_DragDrop(object sender, DragEventArgs e)
        {
            if (e.Data.GetDataPresent(DataFormats.FileDrop))
            {
                string[] files = (string[])e.Data.GetData(DataFormats.FileDrop);
                listBox1.Items.AddRange(files);
            }
        }

        private void browseTargetInput_FileOk(object sender, CancelEventArgs e)
        {

        }

        private void button1_Click(object sender, EventArgs e)
        {
            listBox1.Items.Clear();
        }
    }
}
