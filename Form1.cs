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
using System.Threading;

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
        public static Dictionary<string, List<MSDataScan>> MS1SPECTRA;
        public static Dictionary<int, double> CORRECTIONFACTORS;
        public static double MAXIMUMDNL;
        public static double THEORETICALSEPARATION;
        public static double QUANTRESOLUTION;
        public static double TOLERANCE;
        public static List<NamedChemicalFormula> ISOTOPOLOGUELABELS;
        public static List<NamedChemicalFormula> CLUSTERLABELS;

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
        public static bool CHECKPARTIALINCORPORATION;
        public static bool SEGMENTEDINJECTIONTIMES;
        public static bool OUTPUTRT;
        public static bool LYSINEPURITYCORRECTION;
        public static bool FUSION;
        public static bool OUTPUTSPACINGS;
        public static bool TAGQUANTOUTPUT;
        public static bool OUTPUTPATTERN;
        //public static Dictionary<int, List<double>> PEAKSNPAIRS;
        //public static int QUANTCOUNT;
        //public static bool CALCIUM;
        
        // Experiment types
        //public static bool NEUCODE;
        //public static bool TRADITIONAL;
        //public static bool ICAT;
        //public static bool NEUCODE_DUPLEX_LYS8_36MDA;
        //public static bool NEUCODE_DUPLEX_LEU7_18MDA;
        //public static bool NEUCODE_TRIPLEX_LYS8_18MDA;
        //public static bool NEUCODE_FOURPLEX_LYS8_12MDA;
        //public static bool NEUCODE_SIXPLEX_LYS8_6MDA;
        //public static bool NEUCODE_DUPLEX_LYS1_6MDA;
        //public static bool NEUCODE_DUPLEX_CARBAMYL;
        //public static bool SILAC_DUPLEX_LYSC;
        //public static bool SILAC_DUPLEX_LYSCN;
        //public static bool SILAC_DUPLEX_LYSH;
        //public static bool SILAC_DUPLEX_LEUCN;
        //public static bool SILAC_DUPLEX_LEUH;
        //public static bool NEUCODE_SIXPLEX_MTRAQ;
        //public static bool NEUCODE_SIXPLEX_ARG;
        //public static bool NEUCODE_SIXPLEX_LEU;
        //public static bool NEUCODE_SIXPLEX_DIMETHYL;
        //public static bool NEUCODE_9PLEX_DIMETHYL;
        //public static bool NEUCODE_12PLEX_DIMETHYL;
        //public static bool NEUCODE_18PLEX_DIMETHYL;
        //public static bool NEUCODE_12PLEX;
        //public static bool NEUCODE_4PLEX_HEAVY;
        //public static bool NEUCODE_4PLEX_MEDIUM;
        //public static bool NEUCODE_4PLEX_LIGHT;
        //public static bool NEUCODE_9PLEX_MTRAQ;
        //public static bool NEUCODE_12PLEX_MTRAQ;
        //public static bool NEUCODE_18PLEX_MTRAQ;
        public static int NUMCHANNELS;
        public static int NUMISOTOPOLOGUES;
        public static int NUMCLUSTERS;
        public static bool NHSCLUSTER;
        public static bool LYSCLUSTER;
        public static bool ARGCLUSTER;
        public static bool LEUCLUSTER;
        public static bool CYSCLUSTER;
        public static bool LYSISOTOPOLOGUE;
        public static bool LEUISOTOPOLOGUE;
        public static bool NHSISOTOPOLOGUE;
        public static bool GYGIDIMETHYL;

        public event EventHandler<MessageEventArgs> OnMessage;

        public Form1()
        {

            InitializeComponent();

            rawFileBox.Text = @"E:\Desktop\4-plex Yeast\Optimization\RAW\Duplex Ratios";
            outputFolderBox.Text = @"E:\Desktop\4-plex Yeast\Optimization\QUANT\Duplex Ratios\36 mDa";
            //listBox1.Items.Add(@"E:\Desktop\4-plex Yeast\Optimization\FDR\Individual Incorporation\25Aug2013-Avg\target-decoy\Aug8th2013_ASH_NeuCode_2Plex_30k_480K_DDTop20_12mDa_ITMS_CID_target.csv");
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
            OUTPUTRT = TrackRT.Checked;
            LYSINEPURITYCORRECTION = PurityCorrection.Checked;
            NOISEBANDCAP = noiseBandCap.Checked;
            PEAKCOALESCENCE = coalescence.Checked;
            MULTIINJECT = MultipleInjections.Checked;
            AGCBINNING = AGCBins.Checked;
            FUSION = Fusion.Checked;
            OUTPUTSPACINGS = rawDataOutput.Checked;
            TAGQUANTOUTPUT = tagQuantOutput.Checked;
            NHSCLUSTER = false;
            LYSCLUSTER = false;
            ARGCLUSTER = false;
            LEUCLUSTER = false;
            CYSCLUSTER = false;
            LYSISOTOPOLOGUE = false;
            LEUISOTOPOLOGUE = false;
            NHSISOTOPOLOGUE = false;
            OUTPUTPATTERN = false;

            Dictionary<string, PeptideID> allPeptides = new Dictionary<string, PeptideID>();
            List<PrecursorPPM> PRECURSORPPM = new List<PrecursorPPM>();
            RAWFILES = new Dictionary<string, MSDataFile>();
            List<CoalescenceCheck> INTENSITY_MISSINGCHANNEL = new List<CoalescenceCheck>();
            List<Spacing> spacings = null;
            if (OUTPUTSPACINGS)
            {
                spacings = new List<Spacing>();
            }

            //NEUCODE_DUPLEX_LYS8_36MDA = NeuCodeLys8Duplex.Checked && !mTRAQlight.Checked && !Arg000.Checked && !Leu000.Checked;
            //NEUCODE_DUPLEX_LEU7_18MDA = NeuCodeLeu7Duplex.Checked;
            //NEUCODE_TRIPLEX_LYS8_18MDA = NeuCodeLys8Triplex.Checked && !mTRAQlight.Checked;
            //NEUCODE_FOURPLEX_LYS8_12MDA = NeuCodeLys8Fourplex.Checked && !mTRAQlight.Checked;
            //NEUCODE_SIXPLEX_LYS8_6MDA = NeuCodeLys8Sixplex.Checked && !mTRAQlight.Checked;
            //NEUCODE_9PLEX_MTRAQ = NeuCodeLys8Triplex.Checked && mTRAQlight.Checked;
            //NEUCODE_12PLEX_MTRAQ = NeuCodeLys8Fourplex.Checked && mTRAQlight.Checked;
            //NEUCODE_18PLEX_MTRAQ = NeuCodeLys8Sixplex.Checked && mTRAQlight.Checked;
            //NEUCODE_DUPLEX_LYS1_6MDA = NeuCodeLys1.Checked;
            //NEUCODE_DUPLEX_CARBAMYL = CarbamylCN.Checked;
            //NEUCODE_SIXPLEX_MTRAQ = NeuCodeLys8Duplex.Checked && mTRAQlight.Checked;
            //NEUCODE_SIXPLEX_ARG = NeuCodeLys8Duplex.Checked && Arg000.Checked;
            //NEUCODE_SIXPLEX_LEU = NeuCodeLys8Duplex.Checked && Leu000.Checked;
            //SILAC_DUPLEX_LYSC = SILACLys6C.Checked;
            //SILAC_DUPLEX_LYSCN = SILACLys8CN.Checked;
            //SILAC_DUPLEX_LYSH = SILACLys8D.Checked;
            //SILAC_DUPLEX_LEUCN = SILACLeu7CN.Checked;
            //SILAC_DUPLEX_LEUH = SILACLeu7D.Checked;
            //NEUCODE_12PLEX = Twelveplex.Checked;
            //NEUCODE_4PLEX_HEAVY = FourplexH.Checked;
            //NEUCODE_4PLEX_MEDIUM = FourplexM.Checked;
            //NEUCODE_4PLEX_LIGHT = FourplexL.Checked;
            //ICAT = Icat.Checked;

            setExperimentConfiguration();

            WriteMessage("starting");

            // Cycle through .csv files to make a list of identified peptides and properties
            readCsvInputFile(allPeptides, file);   
      
            int rawFileCount = 0;
            int totalRawFiles = RAWFILES.Count;
            int THEORETICALLYRESOLVABLE = 0;
            int CONTAINSLABEL = 0;

            WriteMessage("calculating systematic error");

            foreach (PeptideID peptide in allPeptides.Values)
            {
                peptide.maximizeResolvability();
                if (peptide.numIsotopologues > 1 && peptide.numIsotopologueLabels > 0)
                {
                    CONTAINSLABEL++;
                    if (peptide.theoreticallyResolvable)
                    {
                        THEORETICALLYRESOLVABLE++;
                    }
                }
                else if (peptide.numIsotopologues < 2 && peptide.numClusterLabels > 0)
                {
                    CONTAINSLABEL++;
                    THEORETICALLYRESOLVABLE++;
                }
            }

            // Calculate systematic error by looking for monoisotopes in MS scan preceding best MS/MS scan
            foreach (MSDataFile rawFile in RAWFILES.Values)
            {
                rawFileCount++;
                rawFile.Open();
                WriteMessage("calculating ppm error in raw file " + rawFileCount + " of " + totalRawFiles);
                
                foreach (PeptideID peptide in allPeptides.Values)
                {
                    if (peptide.theoreticallyResolvable)
                    {
                        PeptideSpectralMatch psm = null;
                        if (peptide.bestPSMs.TryGetValue(rawFile, out psm) && psm.Equals(peptide.bestPSM))
                        {
                            peptide.precursorPPMError(rawFile, PRECURSORPPM);
                        }
                    }
                }
                rawFile.Dispose();
            }
            
            PRECURSORPPM.Sort();
            //writePrecursorPPMOutput(PRECURSORPPM);

            // Only perform systematic error adjustment if 10% or more of peptides are detected
            int peptideVersions = NUMISOTOPOLOGUES;
            if (NUMISOTOPOLOGUES > 2) peptideVersions = NUMCHANNELS;
            else if (NUMCLUSTERS > 1) peptideVersions = NUMISOTOPOLOGUES * NUMCLUSTERS;
            if ((PRECURSORPPM.Count / peptideVersions) > (0.1 * THEORETICALLYRESOLVABLE))
            {
                SYSTEMATICERROR = PRECURSORPPM.ElementAt(PRECURSORPPM.Count / 2).Ppm; //Set systematic error as the median value of precursors
            }
            else
            {
                SYSTEMATICERROR = 0.0;
            }

            WriteMessage("labeled peptides: " + CONTAINSLABEL);
            WriteMessage("theoretically resolvable peptides: " + THEORETICALLYRESOLVABLE);
            WriteMessage("peptides found: " + PRECURSORPPM.Count / peptideVersions);
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
                    if (uniquePeptide.theoreticallyResolvable)
                    {
                        // Only consider peptides that were identified in this raw file
                        List<PeptideSpectralMatch> psms = null;
                        if (uniquePeptide.PSMs.TryGetValue(rawFile, out psms))
                        {
                            uniquePeptide.calculateScanRange(rawFile, RTWINDOW);
                            foreach (MSDataScan currentScan in uniquePeptide.fullScanList)
                            {                               
                                uniquePeptide.findPeaks(currentScan, rawFile, spacings);
                            }

                            List<Pair> pairsFound = uniquePeptide.allPairs[rawFile];
                            if (pairsFound.Count > 0)
                            {
                                uniquePeptide.checkPairSpacing(rawFile);
                                List<Pair> incompletePairsSpacingsChecked = uniquePeptide.allPairs[rawFile];
                                List<Pair> completePairsSpacingsChecked = uniquePeptide.completePairs[rawFile];
                                if (incompletePairsSpacingsChecked.Count + completePairsSpacingsChecked.Count < 1)
                                {
                                    uniquePeptide.noQuantReason = NonQuantifiableType.ImproperSpacing;
                                }
                                else
                                {
                                    if (NOISEBANDCAP)
                                    {
                                        uniquePeptide.sortPairs(rawFile);
                                    }
                                }   
                            }
                            else
                            {
                                uniquePeptide.noQuantReason = NonQuantifiableType.NoPeaksFoundWithinTolerance;
                            }
                        }
                    }
                }
                rawFile.Dispose();
            }

            if (OUTPUTSPACINGS)
            {
                writeSpacingOutput(spacings, file);
            }

            if (NOISEBANDCAP && PEAKCOALESCENCE && NUMISOTOPOLOGUES > 1) // Only functional for 2-plex NeuCode
            {
                WriteMessage("calculating coalescence threshold");
                MAXIMUMDNL = calculateCoalescenceThreshold(file, allPeptides, 0.1, INTENSITY_MISSINGCHANNEL);
                WriteMessage("intensity threshold: " + MAXIMUMDNL);
            }

            // Validate NeuCode pairs with missing channel(s) -- apply noise level or discard due to coalescence
            if (NOISEBANDCAP)
            {
                rawFileCount = 0;
                foreach (MSDataFile rawFile in RAWFILES.Values)
                {
                    rawFileCount++;
                    rawFile.Open();
                    WriteMessage("imputing missing channels in raw file " + rawFileCount + " of " + totalRawFiles);

                    foreach (PeptideID uniquePeptide in allPeptides.Values)
                    {
                        if (uniquePeptide.theoreticallyResolvable)
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
                    rawFile.Dispose();
                }
            }

            WriteMessage("quantifying NeuCode pairs");

            foreach (PeptideID uniquePeptide in allPeptides.Values)
            {
                if (uniquePeptide.theoreticallyResolvable)
                {
                    uniquePeptide.quantify();
                }
                else
                {
                    if (uniquePeptide.numLabels == 0)
                    {
                        uniquePeptide.noQuantReason = NonQuantifiableType.NoLabel;
                    }
                    else
                    {
                        uniquePeptide.noQuantReason = NonQuantifiableType.NotResolvable;
                    }
                    
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
            if (TAGQUANTOUTPUT)
            {
                writeTagQuantOutputFile(allPeptides, file);
            }
            WriteMessage("finished");
        }

        public void WriteMessage(string message)
        {
            if (InvokeRequired)
            {
                Invoke(new OnMessageDelegate(WriteMessage), new object[] { message });
            }
            else
            {
                Console.WriteLine(message);
                richTextBox1.AppendText(message + "\n");
                richTextBox1.ScrollToCaret();                
            }           
        }

        private void writeSpacingOutput(List<Spacing> spacings, string file)
        {
            // Pair spacing outputs can be printed out to a file to assess spacing distributions
            StreamWriter spacingWriter = new StreamWriter(Path.Combine(outputFolderBox.Text, Path.GetFileNameWithoutExtension(file) + "_spacings.csv"));
            string header1 = ("Peptide, RT, Charge, # Labels, Isotope, Light MZ, Light Intensity, Heavy MZ, Heavy Intensity, Theo Spacing (Da), Spacing (Da), Theo Spacing (Th), Spacing (Th)");
            spacingWriter.WriteLine(header1);
            WriteMessage("writing spacing output");

            foreach (Spacing spacing in spacings)
            {
                spacingWriter.WriteLine("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12}", spacing.Sequence, spacing.RetentionTime, spacing.Charge, spacing.NumLabels, spacing.Isotope, spacing.LightMZ, spacing.LightInt, spacing.HeavyMZ, spacing.HeavyInt, spacing.TheoSpacingDa, spacing.SpacingDa, spacing.TheoSpacingTh, spacing.SpacingTh);
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
            // Amino acids
            NamedChemicalFormula K000 = NamedChemicalFormula.AddModification("", "Lys +0");
            NamedChemicalFormula K0001 = NamedChemicalFormula.AddModification("O-1 O{18}1", "Lys +2 18O");
            NamedChemicalFormula K600 = NamedChemicalFormula.AddModification("C-6 C{13}6", "Lys +6 13C6");
            NamedChemicalFormula K602 = NamedChemicalFormula.AddModification("C-6 C{13}6 N-2 N{15}2", "Lys +8 13C6 15N2");
            NamedChemicalFormula K422 = NamedChemicalFormula.AddModification("C-4 C{13}4 H-2 H{2}2 N-2 N{15}2", "Lys +8 13C4 2H2 15N2");
            NamedChemicalFormula K521 = NamedChemicalFormula.AddModification("C-5 C{13}5 H-2 H{2}2 N-1 N{15}1", "Lys +8 13C5 2H2 15N1");
            NamedChemicalFormula K341 = NamedChemicalFormula.AddModification("C-3 C{13}3 H-4 H{2}4 N-1 N{15}1", "Lys +8 13C3 2H4 15N1");
            NamedChemicalFormula K440 = NamedChemicalFormula.AddModification("C-4 C{13}4 H-4 H{2}4", "Lys +8 13C4 2H4");
            NamedChemicalFormula K080 = NamedChemicalFormula.AddModification("H-8 H{2}8", "Lys +8 2H8");
            NamedChemicalFormula R000 = NamedChemicalFormula.AddModification("", "Arg +0");
            NamedChemicalFormula R600 = NamedChemicalFormula.AddModification("C-6 C{13}6", "Arg +6 13C6");
            NamedChemicalFormula R604 = NamedChemicalFormula.AddModification("C-6 C{13}6 N-4 N{15}4", "Arg +10 13C6 15N4");
            NamedChemicalFormula K100 = NamedChemicalFormula.AddModification("C-1 C{13}1", "Lys +1 13C1");
            NamedChemicalFormula K001 = NamedChemicalFormula.AddModification("N-1 N{15}1", "Lys +1 15N1");
            NamedChemicalFormula K040 = NamedChemicalFormula.AddModification("H-4 H{2}4", "Lys +4 2H4");
            NamedChemicalFormula L000 = NamedChemicalFormula.AddModification("", "Leu +0");
            NamedChemicalFormula L601 = NamedChemicalFormula.AddModification("N-1 N{15}1 C-6 C{13}6", "Leu +7 13C6 15N1");
            NamedChemicalFormula L070 = NamedChemicalFormula.AddModification("H-7 H{2}7", "Leu +7 2H7");
            NamedChemicalFormula L600 = NamedChemicalFormula.AddModification("C-6 C{13}6", "Leu +6 13C6");
            NamedChemicalFormula L0100 = NamedChemicalFormula.AddModification("H-10 H{2}10", "Leu +10 2H10");

            // mTRAQ labels
            NamedChemicalFormula lightmTRAQ = NamedChemicalFormula.AddModification("H{1}12 C{12}7 N{14}2 O{16}1", "mTRAQ L");
            NamedChemicalFormula lightmTRAQK602 = NamedChemicalFormula.AddModification("H{1}12 C{12}7 N{14}2 O{16}1 C-6 C{13}6 N-2 N{15}2", "mTRAQ L Lys +8 13C6 15N2");
            NamedChemicalFormula lightmTRAQK422 = NamedChemicalFormula.AddModification("H{1}12 C{12}7 N{14}2 O{16}1 C-4 C{13}4 H-2 H{2}2 N-2 N{15}2", "mTRAQ L Lys +8 13C4 2H2 15N2");
            NamedChemicalFormula lightmTRAQK521 = NamedChemicalFormula.AddModification("H{1}12 C{12}7 N{14}2 O{16}1 C-5 C{13}5 H-2 H{2}2 N-1 N{15}1", "mTRAQ L Lys +8 13C5 2H2 15N1");
            NamedChemicalFormula lightmTRAQK341 = NamedChemicalFormula.AddModification("H{1}12 C{12}7 N{14}2 O{16}1 C-3 C{13}3 H-4 H{2}4 N-1 N{15}1", "mTRAQ L Lys +8 13C3 2H4 15N1");
            NamedChemicalFormula lightmTRAQK440 = NamedChemicalFormula.AddModification("H{1}12 C{12}7 N{14}2 O{16}1 C-4 C{13}4 H-4 H{2}4", "mTRAQ L Lys +8 13C4 2H4");
            NamedChemicalFormula lightmTRAQK080 = NamedChemicalFormula.AddModification("H{1}12 C{12}7 N{14}2 O{16}1 H-8 H{2}8", "mTRAQ L Lys +8 2H8");
            NamedChemicalFormula mediummTRAQ = NamedChemicalFormula.AddModification("H{1}12 C{12}4 C{13}3 N{14}1 N{15}1 O{16}1", "mTRAQ M");
            NamedChemicalFormula mediummTRAQK602 = NamedChemicalFormula.AddModification("H{1}12 C{12}4 C{13}3 N{14}1 N{15}1 O{16}1 C-6 C{13}6 N-2 N{15}2", "mTRAQ M Lys +8 13C6 15N2");
            NamedChemicalFormula mediummTRAQK422 = NamedChemicalFormula.AddModification("H{1}12 C{12}4 C{13}3 N{14}1 N{15}1 O{16}1 C-4 C{13}4 H-2 H{2}2 N-2 N{15}2", "mTRAQ M Lys +8 13C4 2H2 15N2");
            NamedChemicalFormula mediummTRAQK521 = NamedChemicalFormula.AddModification("H{1}12 C{12}4 C{13}3 N{14}1 N{15}1 O{16}1 C-5 C{13}5 H-2 H{2}2 N-1 N{15}1", "mTRAQ M Lys +8 13C5 2H2 15N1");
            NamedChemicalFormula mediummTRAQK341 = NamedChemicalFormula.AddModification("H{1}12 C{12}4 C{13}3 N{14}1 N{15}1 O{16}1 C-3 C{13}3 H-4 H{2}4 N-1 N{15}1", "mTRAQ M Lys +8 13C3 2H4 15N1");
            NamedChemicalFormula mediummTRAQK440 = NamedChemicalFormula.AddModification("H{1}12 C{12}4 C{13}3 N{14}1 N{15}1 O{16}1 C-4 C{13}4 H-4 H{2}4", "mTRAQ M Lys +8 13C4 2H4");
            NamedChemicalFormula mediummTRAQK080 = NamedChemicalFormula.AddModification("H{1}12 C{12}4 C{13}3 N{14}1 N{15}1 O{16}1 H-8 H{2}8", "mTRAQ M Lys +8 2H8");
            NamedChemicalFormula heavymTRAQ = NamedChemicalFormula.AddModification("H{1}12 C{12}1 C{13}6 N{15}2 O{16}1", "mTRAQ H");
            NamedChemicalFormula heavymTRAQK602 = NamedChemicalFormula.AddModification("H{1}12 C{12}1 C{13}6 N{15}2 O{16}1 C-6 C{13}6 N-2 N{15}2", "mTRAQ H Lys +8 13C6 15N2");
            NamedChemicalFormula heavymTRAQK422 = NamedChemicalFormula.AddModification("H{1}12 C{12}1 C{13}6 N{15}2 O{16}1 C-4 C{13}4 H-2 H{2}2 N-2 N{15}2", "mTRAQ H Lys +8 13C4 2H2 15N2");
            NamedChemicalFormula heavymTRAQK521 = NamedChemicalFormula.AddModification("H{1}12 C{12}1 C{13}6 N{15}2 O{16}1 C-5 C{13}5 H-2 H{2}2 N-1 N{15}1", "mTRAQ H Lys +8 13C5 2H2 15N1");
            NamedChemicalFormula heavymTRAQK341 = NamedChemicalFormula.AddModification("H{1}12 C{12}1 C{13}6 N{15}2 O{16}1 C-3 C{13}3 H-4 H{2}4 N-1 N{15}1", "mTRAQ H Lys +8 13C3 2H4 15N1");
            NamedChemicalFormula heavymTRAQK440 = NamedChemicalFormula.AddModification("H{1}12 C{12}1 C{13}6 N{15}2 O{16}1 C-4 C{13}4 H-4 H{2}4", "mTRAQ H Lys +8 13C4 2H4");
            NamedChemicalFormula heavymTRAQK080 = NamedChemicalFormula.AddModification("H{1}12 C{12}1 C{13}6 N{15}2 O{16}1 H-8 H{2}8", "mTRAQ H Lys +8 2H8");

            // Dimethyl labels
            NamedChemicalFormula lightDimethyl = NamedChemicalFormula.AddModification("H{1}4 C{12}2", "Dimethyl L");
            NamedChemicalFormula lightDimethylK602 = NamedChemicalFormula.AddModification("H{1}4 C{12}2 C-6 C{13}6 N-2 N{15}2", "Dimethyl L Lys +8 13C6 15N2");
            NamedChemicalFormula lightDimethylK422 = NamedChemicalFormula.AddModification("H{1}4 C{12}2 C-4 C{13}4 H-2 H{2}2 N-2 N{15}2", "Dimethyl L Lys +8 13C4 2H2 15N2");
            NamedChemicalFormula lightDimethylK521 = NamedChemicalFormula.AddModification("H{1}4 C{12}2 C-5 C{13}5 H-2 H{2}2 N-1 N{15}1", "Dimethyl L Lys +8 13C5 2H2 15N1");
            NamedChemicalFormula lightDimethylK341 = NamedChemicalFormula.AddModification("H{1}4 C{12}2 C-3 C{13}3 H-4 H{2}4 N-1 N{15}1", "Dimethyl L Lys +8 13C3 2H4 15N1");
            NamedChemicalFormula lightDimethylK440 = NamedChemicalFormula.AddModification("H{1}4 C{12}2 C-4 C{13}4 H-4 H{2}4", "Dimethyl L Lys +8 13C4 2H4");
            NamedChemicalFormula lightDimethylK080 = NamedChemicalFormula.AddModification("H{1}4 C{12}2 H-8 H{2}8", "Dimethyl L Lys +8 2H8");
            NamedChemicalFormula mediumDimethyl = NamedChemicalFormula.AddModification("H{2}4 C{12}2", "Dimethyl M");
            NamedChemicalFormula mediumDimethylK602 = NamedChemicalFormula.AddModification("H{2}4 C{12}2 C-6 C{13}6 N-2 N{15}2", "Dimethyl M Lys +8 13C6 15N2");
            NamedChemicalFormula mediumDimethylK422 = NamedChemicalFormula.AddModification("H{2}4 C{12}2 C-4 C{13}4 H-2 H{2}2 N-2 N{15}2", "Dimethyl M Lys +8 13C4 2H2 15N2");
            NamedChemicalFormula mediumDimethylK521 = NamedChemicalFormula.AddModification("H{2}4 C{12}2 C-5 C{13}5 H-2 H{2}2 N-1 N{15}1", "Dimethyl M Lys +8 13C5 2H2 15N1");
            NamedChemicalFormula mediumDimethylK341 = NamedChemicalFormula.AddModification("H{2}4 C{12}2 C-3 C{13}3 H-4 H{2}4 N-1 N{15}1", "Dimethyl M Lys +8 13C3 2H4 15N1");
            NamedChemicalFormula mediumDimethylK440 = NamedChemicalFormula.AddModification("H{2}4 C{12}2 C-4 C{13}4 H-4 H{2}4", "Dimethyl M Lys +8 13C4 2H4");
            NamedChemicalFormula mediumDimethylK080 = NamedChemicalFormula.AddModification("H{2}4 C{12}2 H-8 H{2}8", "Dimethyl M Lys +8 2H8");
            NamedChemicalFormula heavyDimethyl = NamedChemicalFormula.AddModification("H-2 H{2}6 C{13}2", "Dimethyl H");
            NamedChemicalFormula heavyDimethylK602 = NamedChemicalFormula.AddModification("H-2 H{2}6 C{13}2 C-6 C{13}6 N-2 N{15}2", "Dimethyl H Lys +8 13C6 15N2");
            NamedChemicalFormula heavyDimethylK422 = NamedChemicalFormula.AddModification("H-2 H{2}6 C{13}2 C-4 C{13}4 H-2 H{2}2 N-2 N{15}2", "Dimethyl H Lys +8 13C4 2H2 15N2");
            NamedChemicalFormula heavyDimethylK521 = NamedChemicalFormula.AddModification("H-2 H{2}6 C{13}2 C-5 C{13}5 H-2 H{2}2 N-1 N{15}1", "Dimethyl H Lys +8 13C5 2H2 15N1");
            NamedChemicalFormula heavyDimethylK341 = NamedChemicalFormula.AddModification("H-2 H{2}6 C{13}2 C-3 C{13}3 H-4 H{2}4 N-1 N{15}1", "Dimethyl H Lys +8 13C3 2H4 15N1");
            NamedChemicalFormula heavyDimethylK440 = NamedChemicalFormula.AddModification("H-2 H{2}6 C{13}2 C-4 C{13}4 H-4 H{2}4", "Dimethyl H Lys +8 13C4 2H4");
            NamedChemicalFormula heavyDimethylK080 = NamedChemicalFormula.AddModification("H-2 H{2}6 C{13}2 H-8 H{2}8", "Dimethyl H Lys +8 2H8");
            NamedChemicalFormula lightGygiDimethyl = NamedChemicalFormula.AddModification("H{1}4 C{13}2", "Gygi Dimethyl L");
            NamedChemicalFormula heavyGygiDimethyl = NamedChemicalFormula.AddModification("H{1}2 H{2}2 C{12}2", "Gygi Dimethyl H");

            // Chemical labels
            NamedChemicalFormula lightASH1 = NamedChemicalFormula.AddModification("C{12}18 H{1}31 N{14}1 O{16}5 N{15}6", "4plex L1");
            NamedChemicalFormula lightASH2 = NamedChemicalFormula.AddModification("C{12}16 H{1}31 N{14}3 O{16}5 N{15}4 C{13}2", "4plex L2");
            NamedChemicalFormula lightASH3 = NamedChemicalFormula.AddModification("C{12}14 H{1}31 N{14}5 O{16}5 N{15}2 C{13}4", "4plex L3");
            NamedChemicalFormula lightASH4 = NamedChemicalFormula.AddModification("C{12}12 H{1}31 N{14}7 O{16}5 C{13}6", "4plex L4");
            NamedChemicalFormula mediumASH1 = NamedChemicalFormula.AddModification("C{12}14 H{1}31 N{14}1 O{16}5 N{15}6 C{13}4", "4plex M1");
            NamedChemicalFormula mediumASH2 = NamedChemicalFormula.AddModification("C{12}12 H{1}31 N{14}3 O{16}5 N{15}4 C{13}6", "4plex M2");
            NamedChemicalFormula mediumASH3 = NamedChemicalFormula.AddModification("C{12}10 H{1}31 N{14}5 O{16}5 N{15}2 C{13}8", "4plex M3");
            NamedChemicalFormula mediumASH4 = NamedChemicalFormula.AddModification("C{12}8 H{1}31 N{14}7 O{16}5 C{13}10", "4plex M4");
            NamedChemicalFormula heavyASH1 = NamedChemicalFormula.AddModification("C{12}10 H{1}31 N{14}1 O{16}5 N{15}6 C{13}8", "4plex H1");
            NamedChemicalFormula heavyASH2 = NamedChemicalFormula.AddModification("C{12}8 H{1}31 N{14}3 O{16}5 N{15}4 C{13}10", "4plex H2");
            NamedChemicalFormula heavyASH3 = NamedChemicalFormula.AddModification("C{12}6 H{1}31 N{14}5 O{16}5 N{15}2 C{13}12", "4plex H3");
            NamedChemicalFormula heavyASH4 = NamedChemicalFormula.AddModification("C{12}4 H{1}31 N{14}7 O{16}5 C{13}14", "4plex H4");
            NamedChemicalFormula lightCarbamyl = NamedChemicalFormula.AddModification("C{12}1 N{15}1 H{1}2 O{16}1", "Carbamyl L");
            NamedChemicalFormula heavyCarbamyl = NamedChemicalFormula.AddModification("C{13}1 N{14}1 H{1}2 O{16}1", "Carbamyl H");
            NamedChemicalFormula lightICAT = NamedChemicalFormula.AddModification("C{12}10 H{1}15 N{14}3 O{16}2", "iCAT L");
            NamedChemicalFormula heavyICAT = NamedChemicalFormula.AddModification("C{13}9 C{12}1 H{1}15 N{14}3 O{16}2", "iCAT H");

            ISOTOPOLOGUELABELS = new List<NamedChemicalFormula>();
            CLUSTERLABELS = new List<NamedChemicalFormula>();
            
            if (LYSINEPURITYCORRECTION)
            {
                CORRECTIONFACTORS = new Dictionary<int, double>();
                CORRECTIONFACTORS.Add(1, 0.9726);
                CORRECTIONFACTORS.Add(2, 0.7365);
                CORRECTIONFACTORS.Add(3, 0.8887);
                CORRECTIONFACTORS.Add(4, 0.9567);
                CORRECTIONFACTORS.Add(5, 0.8938);
                CORRECTIONFACTORS.Add(6, 0.9618);
            }

            int numLysIsotopologues = 0;
            int numLeuIsotopologues = 0;
            int numIsotopologueTypes = 0;

            // Lys isotopologues
            if (Lys100.Checked)
            {
                ISOTOPOLOGUELABELS.Add(K100);
                numLysIsotopologues++;
            }
            if (Lys001.Checked)
            {
                ISOTOPOLOGUELABELS.Add(K001);
                numLysIsotopologues++;
            }
            if (LysO18.Checked)
            {
                ISOTOPOLOGUELABELS.Add(K0001);
                numLysIsotopologues++;
            }
            if (Lys602.Checked)
            {
                ISOTOPOLOGUELABELS.Add(K602);
                numLysIsotopologues++;
            }
            if (Lys422.Checked)
            {
                ISOTOPOLOGUELABELS.Add(K422);
                numLysIsotopologues++;
            }
            if (Lys521.Checked)
            {
                ISOTOPOLOGUELABELS.Add(K521);
                numLysIsotopologues++;
            }
            if (Lys341.Checked)
            {
                ISOTOPOLOGUELABELS.Add(K341);
                numLysIsotopologues++;
            }
            if (Lys440.Checked)
            {
                ISOTOPOLOGUELABELS.Add(K440);
                numLysIsotopologues++;
            }
            if (Lys080.Checked)
            {
                ISOTOPOLOGUELABELS.Add(K080);
                numLysIsotopologues++;
            }
            if (numLysIsotopologues > 0)
            {
                LYSISOTOPOLOGUE = true;
                numIsotopologueTypes++;
            }

            // Leu isotopologues
            if (Leu601.Checked)
            {
                ISOTOPOLOGUELABELS.Add(L601);
                numLeuIsotopologues++;
            }
            if (Leu070.Checked)
            {
                ISOTOPOLOGUELABELS.Add(L070);
                numLeuIsotopologues++;
            }

            if (numLeuIsotopologues > 0)
            {
                LEUISOTOPOLOGUE = true;
                numIsotopologueTypes++;
            }

            if (numIsotopologueTypes > 1)
            {
                WriteMessage("Experiments with multiple amino acid isotopologue types are not supported");
                ISOTOPOLOGUELABELS.Clear();
                Application.Exit();
            }

            int numLysClusters = 0;
            int numArgClusters = 0;
            int numLeuClusters = 0;
            int numMTRAQClusters = 0;
            int numDimethylClusters = 0;
            int numClusterTypes = 0;

            // Lys clusters
            if (numLysIsotopologues == 1)
            {
                CLUSTERLABELS.Add(ISOTOPOLOGUELABELS[0]);
                ISOTOPOLOGUELABELS.Clear();
                numLysClusters++;
            }
            if (Lys000.Checked)
            {
                numLysClusters++;
            }
            if (Lys040.Checked)
            {
                CLUSTERLABELS.Add(K040);
                numLysClusters++;
            }
            if (numLysClusters > 0)
            {
                LYSCLUSTER = true;
                numClusterTypes++;
            }

            // Arg clusters
            if (Arg000.Checked)
            {
                numArgClusters++;
            }
            if (Arg600.Checked)
            {
                CLUSTERLABELS.Add(R600);
                numArgClusters++;
            }
            if (Arg604.Checked)
            {
                CLUSTERLABELS.Add(R604);
                numArgClusters++;
            }

            if (numArgClusters > 0)
            {
                ARGCLUSTER = true;
                numClusterTypes++;
            }

            // Leu clusters
            if (Leu000.Checked)
            {
                numLeuClusters++;
            }
            if (Leu600.Checked)
            {
                CLUSTERLABELS.Add(L600);
                numLeuClusters++;
            }
            if (Leu0100.Checked)
            {
                CLUSTERLABELS.Add(L0100);
                numLeuClusters++;
            }

            if (numLeuClusters > 0)
            {
                LEUCLUSTER = true;
                numClusterTypes++;
            }

            // mTRAQ clusters
            if (mTRAQlight.Checked)
            {
                CLUSTERLABELS.Add(lightmTRAQ);
                numMTRAQClusters++;
            }
            if (mTRAQmedium.Checked)
            {
                CLUSTERLABELS.Add(mediummTRAQ);
                numMTRAQClusters++;
            }
            if (mTRAQheavy.Checked)
            {
                CLUSTERLABELS.Add(heavymTRAQ);
                numMTRAQClusters++;
            }

            if (numMTRAQClusters > 0)
            {
                NHSCLUSTER = true;
                numClusterTypes++;
            }

            // Dimethyl clusters
            if (dimethyllight.Checked)
            {
                CLUSTERLABELS.Add(lightDimethyl);
                numDimethylClusters++;
            }
            if (dimethylmedium.Checked)
            {
                CLUSTERLABELS.Add(mediumDimethyl);
                numDimethylClusters++;
            }
            if (dimethylheavy.Checked)
            {
                CLUSTERLABELS.Add(heavyDimethyl);
                numDimethylClusters++;
            }

            if (numDimethylClusters > 0)
            {
                NHSCLUSTER = true;
                numClusterTypes++;
            }

            if (numClusterTypes > 1)
            {
                WriteMessage("Experiments with multiple cluster types are not supported");
                CLUSTERLABELS.Clear();
                Application.Exit();
            }

            // Chemical 1 cluster
            if (CarbamylCN.Checked)
            {
                ISOTOPOLOGUELABELS.Add(lightCarbamyl);
                ISOTOPOLOGUELABELS.Add(heavyCarbamyl);
                NHSISOTOPOLOGUE = true;

            }

            else if (GygiDimethyl.Checked)
            {
                ISOTOPOLOGUELABELS.Add(lightGygiDimethyl);
                ISOTOPOLOGUELABELS.Add(heavyGygiDimethyl);
                NHSISOTOPOLOGUE = true;
            }

            else if (FourplexL.Checked)
            {
                ISOTOPOLOGUELABELS.Add(lightASH1);
                ISOTOPOLOGUELABELS.Add(lightASH2);
                ISOTOPOLOGUELABELS.Add(lightASH3);
                ISOTOPOLOGUELABELS.Add(lightASH4);
                NHSCLUSTER = true;
                NHSISOTOPOLOGUE = true;

            }

            else if (FourplexM.Checked)
            {
                ISOTOPOLOGUELABELS.Add(mediumASH1);
                ISOTOPOLOGUELABELS.Add(mediumASH2);
                ISOTOPOLOGUELABELS.Add(mediumASH3);
                ISOTOPOLOGUELABELS.Add(mediumASH4);
                NHSCLUSTER = true;
                NHSISOTOPOLOGUE = true;
            }

            else if (FourplexH.Checked)
            {
                ISOTOPOLOGUELABELS.Add(heavyASH1);
                ISOTOPOLOGUELABELS.Add(heavyASH2);
                ISOTOPOLOGUELABELS.Add(heavyASH3);
                ISOTOPOLOGUELABELS.Add(heavyASH4);
                NHSISOTOPOLOGUE = true;
                NHSCLUSTER = true;
            }

            // Chemical 2 clusters
            if (Icat.Checked)
            {
                CLUSTERLABELS.Add(lightICAT);
                CLUSTERLABELS.Add(heavyICAT);
                CYSCLUSTER = true;

                NUMCLUSTERS = 2;
                NUMISOTOPOLOGUES = 1;
            }

            // Chemical 3 clusters
            else if (Twelveplex.Checked)
            {
                ISOTOPOLOGUELABELS.Add(lightASH1);
                ISOTOPOLOGUELABELS.Add(lightASH2);
                ISOTOPOLOGUELABELS.Add(lightASH3);
                ISOTOPOLOGUELABELS.Add(lightASH4);
                ISOTOPOLOGUELABELS.Add(mediumASH1);
                ISOTOPOLOGUELABELS.Add(mediumASH2);
                ISOTOPOLOGUELABELS.Add(mediumASH3);
                ISOTOPOLOGUELABELS.Add(mediumASH4);
                ISOTOPOLOGUELABELS.Add(heavyASH1);
                ISOTOPOLOGUELABELS.Add(heavyASH2);
                ISOTOPOLOGUELABELS.Add(heavyASH3);
                ISOTOPOLOGUELABELS.Add(heavyASH4);
                NHSISOTOPOLOGUE = true;
                NHSCLUSTER = true;

                NUMCLUSTERS = 3;
                NUMISOTOPOLOGUES = 4;
            }

            else
            {
                if (NHSCLUSTER)
                {
                    if (ISOTOPOLOGUELABELS.Count == 0)
                    {
                        NUMCLUSTERS = CLUSTERLABELS.Count;
                        NUMISOTOPOLOGUES = ISOTOPOLOGUELABELS.Count + 1;
                    }
                    else
                    {
                        NUMISOTOPOLOGUES = ISOTOPOLOGUELABELS.Count;
                        NUMCLUSTERS = CLUSTERLABELS.Count;

                        List<NamedChemicalFormula> NHSISOTOPOLOGUELABELS = new List<NamedChemicalFormula>();

                        if (CLUSTERLABELS.Contains(lightmTRAQ))
                        {
                            if (ISOTOPOLOGUELABELS.Contains(K602)) NHSISOTOPOLOGUELABELS.Add(lightmTRAQK602);
                            if (ISOTOPOLOGUELABELS.Contains(K422)) NHSISOTOPOLOGUELABELS.Add(lightmTRAQK422);
                            if (ISOTOPOLOGUELABELS.Contains(K521)) NHSISOTOPOLOGUELABELS.Add(lightmTRAQK521);
                            if (ISOTOPOLOGUELABELS.Contains(K341)) NHSISOTOPOLOGUELABELS.Add(lightmTRAQK341);
                            if (ISOTOPOLOGUELABELS.Contains(K440)) NHSISOTOPOLOGUELABELS.Add(lightmTRAQK440);
                            if (ISOTOPOLOGUELABELS.Contains(K080)) NHSISOTOPOLOGUELABELS.Add(lightmTRAQK080);
                        }
                        if (CLUSTERLABELS.Contains(mediummTRAQ))
                        {
                            if (ISOTOPOLOGUELABELS.Contains(K602)) NHSISOTOPOLOGUELABELS.Add(mediummTRAQK602);
                            if (ISOTOPOLOGUELABELS.Contains(K422)) NHSISOTOPOLOGUELABELS.Add(mediummTRAQK422);
                            if (ISOTOPOLOGUELABELS.Contains(K521)) NHSISOTOPOLOGUELABELS.Add(mediummTRAQK521);
                            if (ISOTOPOLOGUELABELS.Contains(K341)) NHSISOTOPOLOGUELABELS.Add(mediummTRAQK341);
                            if (ISOTOPOLOGUELABELS.Contains(K440)) NHSISOTOPOLOGUELABELS.Add(mediummTRAQK440);
                            if (ISOTOPOLOGUELABELS.Contains(K080)) NHSISOTOPOLOGUELABELS.Add(mediummTRAQK080);
                        }
                        if (CLUSTERLABELS.Contains(heavymTRAQ))
                        {
                            if (ISOTOPOLOGUELABELS.Contains(K602)) NHSISOTOPOLOGUELABELS.Add(heavymTRAQK602);
                            if (ISOTOPOLOGUELABELS.Contains(K422)) NHSISOTOPOLOGUELABELS.Add(heavymTRAQK422);
                            if (ISOTOPOLOGUELABELS.Contains(K521)) NHSISOTOPOLOGUELABELS.Add(heavymTRAQK521);
                            if (ISOTOPOLOGUELABELS.Contains(K341)) NHSISOTOPOLOGUELABELS.Add(heavymTRAQK341);
                            if (ISOTOPOLOGUELABELS.Contains(K440)) NHSISOTOPOLOGUELABELS.Add(heavymTRAQK440);
                            if (ISOTOPOLOGUELABELS.Contains(K080)) NHSISOTOPOLOGUELABELS.Add(heavymTRAQK080);
                        }
                        if (CLUSTERLABELS.Contains(lightDimethyl))
                        {
                            if (ISOTOPOLOGUELABELS.Contains(K602)) NHSISOTOPOLOGUELABELS.Add(lightDimethylK602);
                            if (ISOTOPOLOGUELABELS.Contains(K422)) NHSISOTOPOLOGUELABELS.Add(lightDimethylK422);
                            if (ISOTOPOLOGUELABELS.Contains(K521)) NHSISOTOPOLOGUELABELS.Add(lightDimethylK521);
                            if (ISOTOPOLOGUELABELS.Contains(K341)) NHSISOTOPOLOGUELABELS.Add(lightDimethylK341);
                            if (ISOTOPOLOGUELABELS.Contains(K440)) NHSISOTOPOLOGUELABELS.Add(lightDimethylK440);
                            if (ISOTOPOLOGUELABELS.Contains(K080)) NHSISOTOPOLOGUELABELS.Add(lightDimethylK080);
                        }
                        if (CLUSTERLABELS.Contains(mediumDimethyl))
                        {
                            if (ISOTOPOLOGUELABELS.Contains(K602)) NHSISOTOPOLOGUELABELS.Add(mediumDimethylK602);
                            if (ISOTOPOLOGUELABELS.Contains(K422)) NHSISOTOPOLOGUELABELS.Add(mediumDimethylK422);
                            if (ISOTOPOLOGUELABELS.Contains(K521)) NHSISOTOPOLOGUELABELS.Add(mediumDimethylK521);
                            if (ISOTOPOLOGUELABELS.Contains(K341)) NHSISOTOPOLOGUELABELS.Add(mediumDimethylK341);
                            if (ISOTOPOLOGUELABELS.Contains(K440)) NHSISOTOPOLOGUELABELS.Add(mediumDimethylK440);
                            if (ISOTOPOLOGUELABELS.Contains(K080)) NHSISOTOPOLOGUELABELS.Add(mediumDimethylK080);
                        }
                        if (CLUSTERLABELS.Contains(heavyDimethyl))
                        {
                            if (ISOTOPOLOGUELABELS.Contains(K602)) NHSISOTOPOLOGUELABELS.Add(heavyDimethylK602);
                            if (ISOTOPOLOGUELABELS.Contains(K422)) NHSISOTOPOLOGUELABELS.Add(heavyDimethylK422);
                            if (ISOTOPOLOGUELABELS.Contains(K521)) NHSISOTOPOLOGUELABELS.Add(heavyDimethylK521);
                            if (ISOTOPOLOGUELABELS.Contains(K341)) NHSISOTOPOLOGUELABELS.Add(heavyDimethylK341);
                            if (ISOTOPOLOGUELABELS.Contains(K440)) NHSISOTOPOLOGUELABELS.Add(heavyDimethylK440);
                            if (ISOTOPOLOGUELABELS.Contains(K080)) NHSISOTOPOLOGUELABELS.Add(heavyDimethylK080);
                        }

                        ISOTOPOLOGUELABELS = NHSISOTOPOLOGUELABELS;

                    }
                }
                else
                {
                    NUMCLUSTERS = CLUSTERLABELS.Count + 1;

                    if (ISOTOPOLOGUELABELS.Count == 0)
                    {
                        NUMISOTOPOLOGUES = 1;
                    }
                    else
                    {
                        NUMISOTOPOLOGUES = ISOTOPOLOGUELABELS.Count;
                    }
                    
                }
            }

            NUMCHANNELS = NUMCLUSTERS * NUMISOTOPOLOGUES;
            
            //// 2-plex ICAT
            //if (ICAT)
            //{
            //    NEUCODE = false;
            //    TRADITIONAL = true;
            //    NUMCHANNELS = 2;
            //    NUMCLUSTERS = 2;
            //    NUMISOTOPOLOGUES = 1;
            //}

            //// 2-plex SILAC
            //if (SILAC_DUPLEX_LYSC || SILAC_DUPLEX_LYSCN || SILAC_DUPLEX_LYSH || SILAC_DUPLEX_LEUCN || SILAC_DUPLEX_LEUH)
            //{
            //    NEUCODE = false;
            //    TRADITIONAL = true;
            //    NUMCHANNELS = 2;
            //    NUMCLUSTERS = 2;
            //    NUMISOTOPOLOGUES = 1;

            //    if (IncompleteIncorporation.Checked)
            //    {
            //        CHECKPARTIALINCORPORATION = true;
            //    }

            //    if (SILAC_DUPLEX_LEUCN && Conversion.Checked)
            //    {
            //        CORRECTLEUNLOSS = true;
            //    }
            //}

            //// 2-plex NeuCode
            //if (NEUCODE_DUPLEX_LYS8_36MDA || NEUCODE_DUPLEX_LYS1_6MDA || NEUCODE_DUPLEX_CARBAMYL)
            //{
            //    NEUCODE = true;
            //    TRADITIONAL = false;
            //    NUMCHANNELS = 2;
            //    NUMCLUSTERS = 1;
            //    NUMISOTOPOLOGUES = 2;
            //}

            //// 3-plex NeuCode
            //if (NEUCODE_TRIPLEX_LYS8_18MDA)
            //{
            //    NEUCODE = true;
            //    TRADITIONAL = false;
            //    NUMCHANNELS = 3;
            //    NUMCLUSTERS = 1;
            //    NUMISOTOPOLOGUES = 3;
            //}

            //// 4-plex NeuCode
            //if (NEUCODE_FOURPLEX_LYS8_12MDA || NEUCODE_4PLEX_LIGHT || NEUCODE_4PLEX_MEDIUM || NEUCODE_4PLEX_HEAVY)
            //{
            //    NEUCODE = true;
            //    TRADITIONAL = false;
            //    NUMCHANNELS = 4;
            //    NUMCLUSTERS = 1;
            //    NUMISOTOPOLOGUES = 4;
            //}

            //// 6-plex NeuCode
            //if (NEUCODE_SIXPLEX_LYS8_6MDA || NEUCODE_SIXPLEX_MTRAQ || NEUCODE_SIXPLEX_ARG || NEUCODE_SIXPLEX_LEU)
            //{
            //    NEUCODE = true;
            //    TRADITIONAL = false;
            //    NUMCHANNELS = 6;

            //    if (NEUCODE_SIXPLEX_MTRAQ || NEUCODE_SIXPLEX_ARG || NEUCODE_SIXPLEX_LEU)
            //    {
            //        NUMCLUSTERS = 3;
            //        NUMISOTOPOLOGUES = 2;

            //        if (NEUCODE_SIXPLEX_ARG && Conversion.Checked)
            //        {
            //            CORRECTARGPROLINE = true;
            //        }

            //        if (NEUCODE_SIXPLEX_LEU && Conversion.Checked)
            //        {
            //            CORRECTLEUDLOSS = true;
            //        }
            //    }
            //    else
            //    {
            //        NUMCLUSTERS = 1;
            //        NUMISOTOPOLOGUES = 6;

            //        if (AGCBINNING)
            //        {
            //            SEGMENTEDINJECTIONTIMES = true;
            //            INJECTIONTIMES = new Dictionary<string, Dictionary<int, Dictionary<Range<double>, double>>>();
            //        }
            //    }
            //}

            //// 9-plex NeuCode
            //if (NEUCODE_9PLEX_MTRAQ)
            //{
            //    NEUCODE = true;
            //    TRADITIONAL = false;
            //    NUMCHANNELS = 9;
            //    NUMCLUSTERS = 3;
            //    NUMISOTOPOLOGUES = 3;
            //}

            //// 12-plex NeuCode
            //if (NEUCODE_12PLEX || NEUCODE_12PLEX_MTRAQ)
            //{
            //    NEUCODE = true;
            //    TRADITIONAL = false;
            //    NUMCHANNELS = 12;
            //    NUMCLUSTERS = 3;
            //    NUMISOTOPOLOGUES = 4;
            //}

            //// 18-plex NeuCode
            //if (NEUCODE_18PLEX_MTRAQ)
            //{
            //    NEUCODE = true;
            //    TRADITIONAL = false;
            //    NUMCHANNELS = 18;
            //    NUMCLUSTERS = 3;
            //    NUMISOTOPOLOGUES = 6;
            //}
        }
        
        private void writeCsvOutputFile(Dictionary<string, PeptideID> allPeptides, string file)
        {
            string outputName = Path.Combine(outputFolderBox.Text, Path.GetFileNameWithoutExtension(file) + "_NeuCode_Quant.csv");
            StreamWriter writer1 = new StreamWriter(outputName);

            if (NUMCHANNELS == 2)
            {
                string header1;
                if (NUMISOTOPOLOGUES > 1)
                {
                    if (OUTPUTRT)
                    {
                        header1 = ("Scan number(s), Raw File(s), Charge State(s), Best Scan Number, Peptide, # Labels, Avg Theo MZ, Avg Adjusted Theo MZ, Resolvable?, Quantified MS1 Scans, # Measurements, Intensity 1, Intensity 2, RT 1, RT 2, Ratio 2/1, Ratio Count, Missing Channels?, Quant Reason");
                        writer1.WriteLine(header1);

                        foreach (PeptideID peptide in allPeptides.Values)
                        {
                            writer1.WriteLine("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18}",
                                peptide.scanNumbers, peptide.rawFiles, peptide.chargeStates, peptide.bestPSM.ScanNumber, peptide.sequence, peptide.numLabels,
                                peptide.averageTheoMZ[0], peptide.averageAdjustedTheoMZ[0],
                                peptide.GetTheoreticalResolvability(0),
                                peptide.countAllPairs, peptide.countAllIsotopes[0], peptide.totalIntensity[0, NUMISOTOPES], peptide.totalIntensity[1, NUMISOTOPES],
                                peptide.maxCompleteIntensity[0, 1], peptide.maxCompleteIntensity[1, 1],
                                peptide.heavyToLightRatioSum[0, 0], peptide.finalQuantified[0], peptide.quantifiedNoiseIncluded[0], peptide.noQuantReason);
                        }
                    }
                    else if (OUTPUTPATTERN)
                    {
                        header1 = ("Scan number(s), Raw File(s), Charge State(s), Best Scan Number, Peptide, # Labels, Avg Theo MZ, Avg Adjusted Theo MZ, Resolvable?, Quantified MS1 Scans, # Measurements, Intensity 1, Intensity 2, Ratio 2/1, Ratio Count, Missing Channels?, Quant Reason, Pattern # Peaks, Pattern ID");
                        writer1.WriteLine(header1);

                        foreach (PeptideID peptide in allPeptides.Values)
                        {
                            writer1.WriteLine("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18}",
                                peptide.scanNumbers, peptide.rawFiles, peptide.chargeStates, peptide.bestPSM.ScanNumber, peptide.sequence, peptide.numLabels,
                                peptide.averageTheoMZ[0], peptide.averageAdjustedTheoMZ[0],
                                peptide.GetTheoreticalResolvability(0),
                                peptide.countAllPairs, peptide.countAllIsotopes[0], peptide.totalIntensity[0, NUMISOTOPES], peptide.totalIntensity[1, NUMISOTOPES],
                                peptide.heavyToLightRatioSum[0, 0], peptide.finalQuantified[0], peptide.quantifiedNoiseIncluded[0], peptide.noQuantReason, peptide.missingChannelPattern[0, 0], peptide.missingChannelPattern[0,1]);
                        }
                    }
                    else
                    {
                        header1 = ("Scan number(s), Raw File(s), Charge State(s), Best Scan Number, Peptide, # Labels, Avg Theo MZ, Avg Adjusted Theo MZ, Resolvable?, Quantified MS1 Scans, # Measurements, Intensity 1, Intensity 2, Ratio 2/1, Ratio Count, Missing Channels?, Quant Reason");
                        writer1.WriteLine(header1);

                        foreach (PeptideID peptide in allPeptides.Values)
                        {
                            writer1.WriteLine("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16}",
                                peptide.scanNumbers, peptide.rawFiles, peptide.chargeStates, peptide.bestPSM.ScanNumber, peptide.sequence, peptide.numLabels,
                                peptide.averageTheoMZ[0], peptide.averageAdjustedTheoMZ[0],
                                peptide.GetTheoreticalResolvability(0),
                                peptide.countAllPairs, peptide.countAllIsotopes[0], peptide.totalIntensity[0, NUMISOTOPES], peptide.totalIntensity[1, NUMISOTOPES],
                                peptide.heavyToLightRatioSum[0, 0], peptide.finalQuantified[0], peptide.quantifiedNoiseIncluded[0], peptide.noQuantReason);
                        }
                    }
                }
                else
                {
                    if (OUTPUTRT)
                    {
                        header1 = ("Scan number(s), Raw File(s), Charge State(s), Best Scan Number, Peptide, # Labels, Cluster 1 Avg Theo MZ, Cluster 2 Avg Theo MZ, Cluster 1 Avg Adjusted Theo MZ, Cluster 2 Avg Adjusted Theo MZ, Resolvable?, Quantified MS1 Scans, # Measurements, Intensity 1, Intensity 2, RT 1, RT 2, Ratio 2/1, Ratio Count, Missing Channels?, Quant Reason");
                        writer1.WriteLine(header1);

                        foreach (PeptideID peptide in allPeptides.Values)
                        {
                            writer1.WriteLine("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20}",
                                peptide.scanNumbers, peptide.rawFiles, peptide.chargeStates, peptide.bestPSM.ScanNumber, peptide.sequence, peptide.numLabels,
                                peptide.averageTheoMZ[0], peptide.averageTheoMZ[1], peptide.averageAdjustedTheoMZ[0], peptide.averageAdjustedTheoMZ[1],
                                peptide.GetTheoreticalResolvability(0),
                                peptide.countAllPairs, peptide.countAllIsotopes[0], peptide.totalIntensity[0, NUMISOTOPES], peptide.totalIntensity[1, NUMISOTOPES],
                                peptide.maxCompleteIntensity[0, 1], peptide.maxCompleteIntensity[1, 1],
                                peptide.heavyToLightRatioSum[0, 0], peptide.finalQuantified[0], peptide.quantifiedNoiseIncluded[0], peptide.noQuantReason);
                        }
                    }
                    else if (OUTPUTPATTERN)
                    {
                        header1 = ("Scan number(s), Raw File(s), Charge State(s), Best Scan Number, Peptide, # Labels, Cluster 1 Avg Theo MZ, Cluster 2 Avg Theo MZ, Cluster 1 Avg Adjusted Theo MZ, Cluster 2 Avg Adjusted Theo MZ, Resolvable?, Quantified MS1 Scans, # Measurements, Intensity 1, Intensity 2, Ratio 2/1, Ratio Count, Missing Channels?, Quant Reason, Pattern # Peaks, Pattern ID");
                        writer1.WriteLine(header1);

                        foreach (PeptideID peptide in allPeptides.Values)
                        {
                            writer1.WriteLine("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20}",
                                peptide.scanNumbers, peptide.rawFiles, peptide.chargeStates, peptide.bestPSM.ScanNumber, peptide.sequence, peptide.numLabels,
                                peptide.averageTheoMZ[0], peptide.averageTheoMZ[1], peptide.averageAdjustedTheoMZ[0], peptide.averageAdjustedTheoMZ[1],
                                peptide.GetTheoreticalResolvability(0),
                                peptide.countAllPairs, peptide.countAllIsotopes[0], peptide.totalIntensity[0, NUMISOTOPES], peptide.totalIntensity[1, NUMISOTOPES],
                                peptide.heavyToLightRatioSum[0, 0], peptide.finalQuantified[0], peptide.quantifiedNoiseIncluded[0], peptide.noQuantReason, peptide.missingChannelPattern[0,0], peptide.missingChannelPattern[0,1]);
                        }
                    }
                    else
                    {
                        header1 = ("Scan number(s), Raw File(s), Charge State(s), Best Scan Number, Peptide, # Labels, Cluster 1 Avg Theo MZ, Cluster 2 Avg Theo MZ, Cluster 1 Avg Adjusted Theo MZ, Cluster 2 Avg Adjusted Theo MZ, Resolvable?, Quantified MS1 Scans, # Measurements, Intensity 1, Intensity 2, Ratio 2/1, Ratio Count, Missing Channels?, Quant Reason");
                        writer1.WriteLine(header1);

                        foreach (PeptideID peptide in allPeptides.Values)
                        {
                            writer1.WriteLine("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18}",
                                peptide.scanNumbers, peptide.rawFiles, peptide.chargeStates, peptide.bestPSM.ScanNumber, peptide.sequence, peptide.numLabels,
                                peptide.averageTheoMZ[0], peptide.averageTheoMZ[1], peptide.averageAdjustedTheoMZ[0], peptide.averageAdjustedTheoMZ[1],
                                peptide.GetTheoreticalResolvability(0),
                                peptide.countAllPairs, peptide.countAllIsotopes[0], peptide.totalIntensity[0, NUMISOTOPES], peptide.totalIntensity[1, NUMISOTOPES],
                                peptide.heavyToLightRatioSum[0, 0], peptide.finalQuantified[0], peptide.quantifiedNoiseIncluded[0], peptide.noQuantReason);
                        }
                    }
                }
                writer1.Close();
            }
            else if (NUMCHANNELS == 3)
            {
                string header1;
                if (NUMISOTOPOLOGUES > 1)
                {
                    if (OUTPUTRT)
                    {
                        header1 = ("Scan number(s), Raw File(s), Charge State(s), Best Scan Number, Peptide, # Labels, Avg Theo MZ, Avg Adjusted Theo MZ, Resolvable?, Quantified MS1 Scans, # Measurements, Intensity 1, Intensity 2, Intensity 3, RT 1, RT 2, RT 3, Ratio 2/1, Ratio 3/1, Ratio Count, Missing Channels?, Quant Reason");
                        writer1.WriteLine(header1);

                        foreach (PeptideID peptide in allPeptides.Values)
                        {
                            writer1.WriteLine("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21}",
                                peptide.scanNumbers, peptide.rawFiles, peptide.chargeStates, peptide.bestPSM.ScanNumber, peptide.sequence, peptide.numLabels,
                                peptide.averageTheoMZ[0], peptide.averageAdjustedTheoMZ[0],
                                peptide.GetTheoreticalResolvability(0),
                                peptide.countAllPairs, peptide.countAllIsotopes[0], peptide.totalIntensity[0, NUMISOTOPES], peptide.totalIntensity[1, NUMISOTOPES], peptide.totalIntensity[2, NUMISOTOPES],
                                peptide.maxCompleteIntensity[0, 1], peptide.maxCompleteIntensity[1, 1], peptide.maxCompleteIntensity[2, 1],
                                peptide.heavyToLightRatioSum[0, 0], peptide.heavyToLightRatioSum[1, 0], peptide.finalQuantified[0], peptide.quantifiedNoiseIncluded[0], peptide.noQuantReason);
                        }
                    }
                    else if (OUTPUTPATTERN)
                    {
                        header1 = ("Scan number(s), Raw File(s), Charge State(s), Best Scan Number, Peptide, # Labels, Avg Theo MZ, Avg Adjusted Theo MZ, Resolvable?, Quantified MS1 Scans, # Measurements, Intensity 1, Intensity 2, Intensity 3, Ratio 2/1, Ratio 3/1, Ratio Count, Missing Channels?, Quant Reason, Pattern # Peaks, Pattern ID");
                        writer1.WriteLine(header1);

                        foreach (PeptideID peptide in allPeptides.Values)
                        {
                            writer1.WriteLine("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20}",
                                peptide.scanNumbers, peptide.rawFiles, peptide.chargeStates, peptide.bestPSM.ScanNumber, peptide.sequence, peptide.numLabels,
                                peptide.averageTheoMZ[0], peptide.averageAdjustedTheoMZ[0],
                                peptide.GetTheoreticalResolvability(0),
                                peptide.countAllPairs, peptide.countAllIsotopes[0], peptide.totalIntensity[0, NUMISOTOPES], peptide.totalIntensity[1, NUMISOTOPES], peptide.totalIntensity[2, NUMISOTOPES],
                                peptide.heavyToLightRatioSum[0, 0], peptide.heavyToLightRatioSum[1, 0], peptide.finalQuantified[0], peptide.quantifiedNoiseIncluded[0], peptide.noQuantReason, peptide.missingChannelPattern[0,0], peptide.missingChannelPattern[0,1]);
                        }
                    }
                    else
                    {
                        header1 = ("Scan number(s), Raw File(s), Charge State(s), Best Scan Number, Peptide, # Labels, Avg Theo MZ, Avg Adjusted Theo MZ, Resolvable?, Quantified MS1 Scans, # Measurements, Intensity 1, Intensity 2, Intensity 3, Ratio 2/1, Ratio 3/1, Ratio Count, Missing Channels?, Quant Reason");
                        writer1.WriteLine(header1);

                        foreach (PeptideID peptide in allPeptides.Values)
                        {
                            writer1.WriteLine("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18}",
                                peptide.scanNumbers, peptide.rawFiles, peptide.chargeStates, peptide.bestPSM.ScanNumber, peptide.sequence, peptide.numLabels,
                                peptide.averageTheoMZ[0], peptide.averageAdjustedTheoMZ[0],
                                peptide.GetTheoreticalResolvability(0),
                                peptide.countAllPairs, peptide.countAllIsotopes[0], peptide.totalIntensity[0, NUMISOTOPES], peptide.totalIntensity[1, NUMISOTOPES], peptide.totalIntensity[2, NUMISOTOPES],
                                peptide.heavyToLightRatioSum[0, 0], peptide.heavyToLightRatioSum[1, 0], peptide.finalQuantified[0], peptide.quantifiedNoiseIncluded[0], peptide.noQuantReason);
                        }
                    }
                }
                else
                {
                    if (OUTPUTRT)
                    {
                        header1 = ("Scan number(s), Raw File(s), Charge State(s), Best Scan Number, Peptide, # Labels, Cluster 1 Avg Theo MZ, Cluster 2 Avg Theo MZ, Cluster 3 Avg Theo MZ, Cluster 1 Avg Adjusted Theo MZ, Cluster 2 Avg Adjusted Theo MZ, Cluster 3 Avg Adjusted Theo MZ, Resolvable?, Quantified MS1 Scans, # Measurements, Intensity 1, Intensity 2, Intensity 3, RT 1, RT 2, RT 3, Ratio 2/1, Ratio 3/1, Ratio Count, Missing Channels?, Quant Reason");
                        writer1.WriteLine(header1);

                        foreach (PeptideID peptide in allPeptides.Values)
                        {
                            writer1.WriteLine("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23},{24},{25}",
                                peptide.scanNumbers, peptide.rawFiles, peptide.chargeStates, peptide.bestPSM.ScanNumber, peptide.sequence, peptide.numLabels,
                                peptide.averageTheoMZ[0], peptide.averageTheoMZ[1], peptide.averageTheoMZ[2], peptide.averageAdjustedTheoMZ[0], peptide.averageAdjustedTheoMZ[1], peptide.averageAdjustedTheoMZ[2],
                                peptide.GetTheoreticalResolvability(0),
                                peptide.countAllPairs, peptide.countAllIsotopes[0], peptide.totalIntensity[0, NUMISOTOPES], peptide.totalIntensity[1, NUMISOTOPES], peptide.totalIntensity[2, NUMISOTOPES],
                                peptide.maxCompleteIntensity[0, 1], peptide.maxCompleteIntensity[1, 1], peptide.maxCompleteIntensity[2, 1],
                                peptide.heavyToLightRatioSum[0, 0], peptide.heavyToLightRatioSum[1, 0], peptide.finalQuantified[0], peptide.quantifiedNoiseIncluded[0], peptide.noQuantReason);
                        }
                    }
                    else
                    {
                        header1 = ("Scan number(s), Raw File(s), Charge State(s), Best Scan Number, Peptide, # Labels, Cluster 1 Avg Theo MZ, Cluster 2 Avg Theo MZ, Cluster 3 Avg Theo MZ, Cluster 1 Avg Adjusted Theo MZ, Cluster 2 Avg Adjusted Theo MZ, Cluster 3 Avg Adjusted Theo MZ, Resolvable?, Quantified MS1 Scans, # Measurements, Intensity 1, Intensity 2, Intensity 3, Ratio 2/1, Ratio 3/1, Ratio Count, Missing Channels?, Quant Reason");
                        writer1.WriteLine(header1);

                        foreach (PeptideID peptide in allPeptides.Values)
                        {
                            writer1.WriteLine("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22}",
                                peptide.scanNumbers, peptide.rawFiles, peptide.chargeStates, peptide.bestPSM.ScanNumber, peptide.sequence, peptide.numLabels,
                                peptide.averageTheoMZ[0], peptide.averageTheoMZ[1], peptide.averageTheoMZ[2], peptide.averageAdjustedTheoMZ[0], peptide.averageAdjustedTheoMZ[1], peptide.averageAdjustedTheoMZ[2],
                                peptide.GetTheoreticalResolvability(0),
                                peptide.countAllPairs, peptide.countAllIsotopes[0], peptide.totalIntensity[0, NUMISOTOPES], peptide.totalIntensity[1, NUMISOTOPES], peptide.totalIntensity[2, NUMISOTOPES],
                                peptide.heavyToLightRatioSum[0, 0], peptide.heavyToLightRatioSum[1, 0], peptide.finalQuantified[0], peptide.quantifiedNoiseIncluded[0], peptide.noQuantReason);
                        }
                    }
                }                
                writer1.Close();
            }
            else if (NUMCHANNELS == 4)
            {
                string header1;
                if (NUMCLUSTERS == 1)
                {
                    if (OUTPUTRT)
                    {
                        header1 = ("Scan number(s), Raw File(s), Charge State(s), Best Scan Number, Peptide, # Labels, Avg Theo MZ, Avg Adjusted Theo MZ, Resolvable?, Quantified MS1 Scans, # Measurements, Intensity 1, Intensity 2, Intensity 3, Intensity 4, RT 1, RT 2, RT 3, RT 4, Ratio 2/1, Ratio 3/1, Ratio 4/1, Ratio Count, Missing Channels?, Quant Reason");
                        writer1.WriteLine(header1);
                        foreach (PeptideID peptide in allPeptides.Values)
                        {
                            writer1.WriteLine("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23},{24}",
                                peptide.scanNumbers, peptide.rawFiles, peptide.chargeStates, peptide.bestPSM.ScanNumber, peptide.sequence, peptide.numLabels,
                                peptide.averageTheoMZ[0], peptide.averageAdjustedTheoMZ[0],
                                peptide.GetTheoreticalResolvability(0),
                                peptide.countAllPairs, peptide.countAllIsotopes[0], peptide.totalIntensity[0, NUMISOTOPES], peptide.totalIntensity[1, NUMISOTOPES], peptide.totalIntensity[2, NUMISOTOPES], peptide.totalIntensity[3, NUMISOTOPES],
                                peptide.maxCompleteIntensity[0, 1], peptide.maxCompleteIntensity[1, 1], peptide.maxCompleteIntensity[2, 1], peptide.maxCompleteIntensity[3, 1],
                                peptide.heavyToLightRatioSum[0, 0], peptide.heavyToLightRatioSum[1, 0], peptide.heavyToLightRatioSum[2, 0], peptide.finalQuantified[0], peptide.quantifiedNoiseIncluded[0], peptide.noQuantReason);
                        }
                    }
                    else if (OUTPUTPATTERN)
                    {
                        header1 = ("Scan number(s), Raw File(s), Charge State(s), Best Scan Number, Peptide, # Labels, Avg Theo MZ, Avg Adjusted Theo MZ, Resolvable?, Quantified MS1 Scans, # Measurements, Intensity 1, Intensity 2, Intensity 3, Intensity 4, Ratio 2/1, Ratio 3/1, Ratio 4/1, Ratio Count, Missing Channels?, Quant Reason, Pattern # Peaks, Pattern ID");
                        writer1.WriteLine(header1);
                        foreach (PeptideID peptide in allPeptides.Values)
                        {
                            writer1.WriteLine("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22}",
                                peptide.scanNumbers, peptide.rawFiles, peptide.chargeStates, peptide.bestPSM.ScanNumber, peptide.sequence, peptide.numLabels,
                                peptide.averageTheoMZ[0], peptide.averageAdjustedTheoMZ[0],
                                peptide.GetTheoreticalResolvability(0),
                                peptide.countAllPairs, peptide.countAllIsotopes[0], peptide.totalIntensity[0, NUMISOTOPES], peptide.totalIntensity[1, NUMISOTOPES], peptide.totalIntensity[2, NUMISOTOPES], peptide.totalIntensity[3, NUMISOTOPES],
                                peptide.heavyToLightRatioSum[0, 0], peptide.heavyToLightRatioSum[1, 0], peptide.heavyToLightRatioSum[2, 0], peptide.finalQuantified[0], peptide.quantifiedNoiseIncluded[0], peptide.noQuantReason, peptide.missingChannelPattern[0,0], peptide.missingChannelPattern[0,1]);
                        }
                    }
                    else
                    {
                        header1 = ("Scan number(s), Raw File(s), Charge State(s), Best Scan Number, Peptide, # Labels, Avg Theo MZ, Avg Adjusted Theo MZ, Resolvable?, Quantified MS1 Scans, # Measurements, Intensity 1, Intensity 2, Intensity 3, Intensity 4, Ratio 2/1, Ratio 3/1, Ratio 4/1, Ratio Count, Missing Channels?, Quant Reason");
                        writer1.WriteLine(header1);
                        foreach (PeptideID peptide in allPeptides.Values)
                        {
                            writer1.WriteLine("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20}",
                                peptide.scanNumbers, peptide.rawFiles, peptide.chargeStates, peptide.bestPSM.ScanNumber, peptide.sequence, peptide.numLabels,
                                peptide.averageTheoMZ[0], peptide.averageAdjustedTheoMZ[0],
                                peptide.GetTheoreticalResolvability(0),
                                peptide.countAllPairs, peptide.countAllIsotopes[0], peptide.totalIntensity[0, NUMISOTOPES], peptide.totalIntensity[1, NUMISOTOPES], peptide.totalIntensity[2, NUMISOTOPES], peptide.totalIntensity[3, NUMISOTOPES],
                                peptide.heavyToLightRatioSum[0, 0], peptide.heavyToLightRatioSum[1, 0], peptide.heavyToLightRatioSum[2, 0], peptide.finalQuantified[0], peptide.quantifiedNoiseIncluded[0], peptide.noQuantReason);
                        }
                    }
                }
                else if (NUMCLUSTERS == 2)
                {
                    if (OUTPUTRT)
                    {
                        header1 = ("Scan number(s), Raw File(s), Charge State(s), Best Scan Number, Peptide, # Labels, Cluster 1 Avg Theo MZ, Cluster 2 Avg Theo MZ, Cluster 1 Avg Adjusted Theo MZ, Cluster 2 Avg Adjusted Theo MZ, Cluster 1 Resolvable?, Cluster 2 Resolvable?, Quantified MS1 Scans, Cluster 1 # Measurements, Cluster 2 # Measurements, Intensity 1, Intensity 2, Intensity 3, Intensity 4, RT 1, RT 2, RT 3, RT 4, Ratio 2/1, Ratio 4/3, Cluster 1 Ratio Count, Cluster 2 Ratio Count, Cluster 1 Missing Channels?, Cluster 2 Missing Channels?, Quant Reason");
                        writer1.WriteLine(header1);
                        foreach (PeptideID peptide in allPeptides.Values)
                        {
                            writer1.WriteLine("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23},{24},{25},{26},{27},{28},{29},{30}",
                                peptide.scanNumbers, peptide.rawFiles, peptide.chargeStates, peptide.bestPSM.ScanNumber, peptide.sequence, peptide.numLabels,
                                peptide.averageTheoMZ[0], peptide.averageTheoMZ[1], peptide.averageAdjustedTheoMZ[0], peptide.averageAdjustedTheoMZ[1],
                                peptide.GetTheoreticalResolvability(0), peptide.GetTheoreticalResolvability(1),
                                peptide.countAllPairs, peptide.countAllIsotopes[0], peptide.countAllIsotopes[1], peptide.totalIntensity[0, NUMISOTOPES], peptide.totalIntensity[1, NUMISOTOPES], peptide.totalIntensity[2, NUMISOTOPES], peptide.totalIntensity[3, NUMISOTOPES],
                                peptide.maxCompleteIntensity[0, 1], peptide.maxCompleteIntensity[1, 1], peptide.maxCompleteIntensity[2, 1], peptide.maxCompleteIntensity[3, 1],
                                peptide.heavyToLightRatioSum[0, 0], peptide.heavyToLightRatioSum[1, 0], peptide.finalQuantified[0], peptide.finalQuantified[1], peptide.quantifiedNoiseIncluded[0], peptide.quantifiedNoiseIncluded[1], peptide.noQuantReason);
                        }
                    }
                    else
                    {
                        header1 = ("Scan number(s), Raw File(s), Charge State(s), Best Scan Number, Peptide, # Labels, Cluster 1 Avg Theo MZ, Cluster 2 Avg Theo MZ, Cluster 1 Avg Adjusted Theo MZ, Cluster 2 Avg Adjusted Theo MZ, Cluster 1 Resolvable?, Cluster 2 Resolvable?, Quantified MS1 Scans, Cluster 1 # Measurements, Cluster 2 # Measurements, Intensity 1, Intensity 2, Intensity 3, Intensity 4, Ratio 2/1, Ratio 4/3, Cluster 1 Ratio Count, Cluster 2 Ratio Count, Cluster 1 Missing Channels?, Cluster 2 Missing Channels?, Quant Reason");
                        writer1.WriteLine(header1);
                        foreach (PeptideID peptide in allPeptides.Values)
                        {
                            writer1.WriteLine("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23},{24},{25},{26}",
                                peptide.scanNumbers, peptide.rawFiles, peptide.chargeStates, peptide.bestPSM.ScanNumber, peptide.sequence, peptide.numLabels,
                                peptide.averageTheoMZ[0], peptide.averageTheoMZ[1], peptide.averageAdjustedTheoMZ[0], peptide.averageAdjustedTheoMZ[1],
                                peptide.GetTheoreticalResolvability(0), peptide.GetTheoreticalResolvability(1),
                                peptide.countAllPairs, peptide.countAllIsotopes[0], peptide.countAllIsotopes[1], peptide.totalIntensity[0, NUMISOTOPES], peptide.totalIntensity[1, NUMISOTOPES], peptide.totalIntensity[2, NUMISOTOPES], peptide.totalIntensity[3, NUMISOTOPES],
                                peptide.heavyToLightRatioSum[0, 0], peptide.heavyToLightRatioSum[1, 0], peptide.finalQuantified[0], peptide.finalQuantified[1], peptide.quantifiedNoiseIncluded[0], peptide.quantifiedNoiseIncluded[1], peptide.noQuantReason);
                        }
                    }
                }
                writer1.Close();
            }
            else if (NUMCHANNELS == 6)
            {
                string header1;
                if (NUMCLUSTERS == 1)
                {
                    if (OUTPUTRT)
                    {
                        header1 = ("Scan number(s), Raw File(s), Charge State(s), Best Scan Number, Peptide, # Labels, Avg Theo MZ, Avg Adjusted Theo MZ, Resolvable?, Quantified MS1 Scans, # Measurements, Intensity 1, Intensity 2, Intensity 3, Intensity 4, Intensity 5, Intensity 6, RT 1, RT 2, RT 3, RT 4, RT 5, RT 6, Ratio 2/1, Ratio 3/1, Ratio 4/1, Ratio 5/1, Ratio 6/1, Ratio Count, Missing Channels?, Quant Reason");
                        writer1.WriteLine(header1);
                        foreach (PeptideID peptide in allPeptides.Values)
                        {
                            writer1.WriteLine("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23},{24},{25},{26},{27},{28},{29},{30}",
                                peptide.scanNumbers, peptide.rawFiles, peptide.chargeStates, peptide.bestPSM.ScanNumber, peptide.sequence, peptide.numLabels,
                                peptide.averageTheoMZ[0], peptide.averageAdjustedTheoMZ[0],
                                peptide.GetTheoreticalResolvability(0),
                                peptide.countAllPairs, peptide.countAllIsotopes[0], peptide.totalIntensity[0, NUMISOTOPES], peptide.totalIntensity[1, NUMISOTOPES], peptide.totalIntensity[2, NUMISOTOPES], peptide.totalIntensity[3, NUMISOTOPES], peptide.totalIntensity[4, NUMISOTOPES], peptide.totalIntensity[5, NUMISOTOPES],
                                peptide.maxCompleteIntensity[0, 1], peptide.maxCompleteIntensity[1, 1], peptide.maxCompleteIntensity[2, 1], peptide.maxCompleteIntensity[3, 1], peptide.maxCompleteIntensity[4,1], peptide.maxCompleteIntensity[5, 1],
                                peptide.heavyToLightRatioSum[0, 0], peptide.heavyToLightRatioSum[1, 0], peptide.heavyToLightRatioSum[2, 0], peptide.heavyToLightRatioSum[3,0], peptide.heavyToLightRatioSum[4,0], peptide.finalQuantified[0], peptide.quantifiedNoiseIncluded[0], peptide.noQuantReason);
                        }
                    }
                    else
                    {
                        header1 = ("Scan number(s), Raw File(s), Charge State(s), Best Scan Number, Peptide, # Labels, Avg Theo MZ, Avg Adjusted Theo MZ, Resolvable?, Quantified MS1 Scans, # Measurements, Intensity 1, Intensity 2, Intensity 3, Intensity 4, Intensity 5, Intensity 6, Ratio 2/1, Ratio 3/1, Ratio 4/1, Ratio 5/1, Ratio 6/1, Ratio Count, Missing Channels?, Quant Reason");
                        writer1.WriteLine(header1);
                        foreach (PeptideID peptide in allPeptides.Values)
                        {
                            writer1.WriteLine("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23},{24}",
                                peptide.scanNumbers, peptide.rawFiles, peptide.chargeStates, peptide.bestPSM.ScanNumber, peptide.sequence, peptide.numLabels,
                                peptide.averageTheoMZ[0], peptide.averageAdjustedTheoMZ[0],
                                peptide.GetTheoreticalResolvability(0),
                                peptide.countAllPairs, peptide.countAllIsotopes[0], peptide.totalIntensity[0, NUMISOTOPES], peptide.totalIntensity[1, NUMISOTOPES], peptide.totalIntensity[2, NUMISOTOPES], peptide.totalIntensity[3, NUMISOTOPES], peptide.totalIntensity[4, NUMISOTOPES], peptide.totalIntensity[5, NUMISOTOPES],
                                peptide.heavyToLightRatioSum[0, 0], peptide.heavyToLightRatioSum[1, 0], peptide.heavyToLightRatioSum[2, 0], peptide.heavyToLightRatioSum[3, 0], peptide.heavyToLightRatioSum[4, 0], peptide.finalQuantified[0], peptide.quantifiedNoiseIncluded[0], peptide.noQuantReason);
                        }
                    }
                }
                else if (NUMCLUSTERS == 2)
                {
                    if (OUTPUTRT)
                    {
                        header1 = ("Scan number(s), Raw File(s), Charge State(s), Best Scan Number, Peptide, # Labels, Cluster 1 Avg Theo MZ, Cluster 2 Avg Theo MZ, Cluster 1 Avg Adjusted Theo MZ, Cluster 2 Avg Adjusted Theo MZ, Cluster 1 Resolvable?, Cluster 2 Resolvable?, Quantified MS1 Scans, Cluster 1 # Measurements, Cluster 2 # Measurements, Intensity 1, Intensity 2, Intensity 3, Intensity 4, Intensity 5, Intensity 6, RT 1, RT 2, RT 3, RT 4, RT 5, RT 6, Ratio 2/1, Ratio 3/1, Ratio 5/4, Ratio 6/4, Cluster 1 Ratio Count, Cluster 2 Ratio Count, Cluster 1 Missing Channels?, Cluster 2 Missing Channels?, Quant Reason");
                        writer1.WriteLine(header1);
                        foreach (PeptideID peptide in allPeptides.Values)
                        {
                            writer1.WriteLine("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23},{24},{25},{26},{27},{28},{29},{30},{31},{32},{33},{34},{35}",
                                peptide.scanNumbers, peptide.rawFiles, peptide.chargeStates, peptide.bestPSM.ScanNumber, peptide.sequence, peptide.numLabels,
                                peptide.averageTheoMZ[0], peptide.averageTheoMZ[1], peptide.averageAdjustedTheoMZ[0], peptide.averageAdjustedTheoMZ[1],
                                peptide.GetTheoreticalResolvability(0), peptide.GetTheoreticalResolvability(1),
                                peptide.countAllPairs, peptide.countAllIsotopes[0], peptide.countAllIsotopes[1], peptide.totalIntensity[0, NUMISOTOPES], peptide.totalIntensity[1, NUMISOTOPES], peptide.totalIntensity[2, NUMISOTOPES], peptide.totalIntensity[3, NUMISOTOPES], peptide.totalIntensity[4, NUMISOTOPES], peptide.totalIntensity[5, NUMISOTOPES],
                                peptide.maxCompleteIntensity[0, 1], peptide.maxCompleteIntensity[1, 1], peptide.maxCompleteIntensity[2, 1], peptide.maxCompleteIntensity[3, 1], peptide.maxCompleteIntensity[4, 1], peptide.maxCompleteIntensity[5, 1],
                                peptide.heavyToLightRatioSum[0, 0], peptide.heavyToLightRatioSum[1, 0], peptide.heavyToLightRatioSum[2, 0], peptide.heavyToLightRatioSum[3, 0], peptide.finalQuantified[0], peptide.finalQuantified[1], peptide.quantifiedNoiseIncluded[0], peptide.quantifiedNoiseIncluded[1], peptide.noQuantReason);
                        }
                    }
                    else
                    {
                        header1 = ("Scan number(s), Raw File(s), Charge State(s), Best Scan Number, Peptide, # Labels, Cluster 1 Avg Theo MZ, Cluster 2 Avg Theo MZ, Cluster 1 Avg Adjusted Theo MZ, Cluster 2 Avg Adjusted Theo MZ, Cluster 1 Resolvable?, Cluster 2 Resolvable?, Quantified MS1 Scans, Cluster 1 # Measurements, Cluster 2 # Measurements, Intensity 1, Intensity 2, Intensity 3, Intensity 4, Intensity 5, Intensity 6, Ratio 2/1, Ratio 3/1, Ratio 5/4, Ratio 6/4, Cluster 1 Ratio Count, Cluster 2 Ratio Count, Cluster 1 Missing Channels?, Cluster 2 Missing Channels?, Quant Reason");
                        writer1.WriteLine(header1);
                        foreach (PeptideID peptide in allPeptides.Values)
                        {
                            writer1.WriteLine("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23},{24},{25},{26},{27},{28},{29}",
                                peptide.scanNumbers, peptide.rawFiles, peptide.chargeStates, peptide.bestPSM.ScanNumber, peptide.sequence, peptide.numLabels,
                                peptide.averageTheoMZ[0], peptide.averageTheoMZ[1], peptide.averageAdjustedTheoMZ[0], peptide.averageAdjustedTheoMZ[1],
                                peptide.GetTheoreticalResolvability(0), peptide.GetTheoreticalResolvability(1),
                                peptide.countAllPairs, peptide.countAllIsotopes[0], peptide.countAllIsotopes[1], peptide.totalIntensity[0, NUMISOTOPES], peptide.totalIntensity[1, NUMISOTOPES], peptide.totalIntensity[2, NUMISOTOPES], peptide.totalIntensity[3, NUMISOTOPES], peptide.totalIntensity[4, NUMISOTOPES], peptide.totalIntensity[5, NUMISOTOPES],
                                peptide.heavyToLightRatioSum[0, 0], peptide.heavyToLightRatioSum[1, 0], peptide.heavyToLightRatioSum[2, 0], peptide.heavyToLightRatioSum[3, 0], peptide.finalQuantified[0], peptide.finalQuantified[1], peptide.quantifiedNoiseIncluded[0], peptide.quantifiedNoiseIncluded[1], peptide.noQuantReason);
                        }
                    }
                }
                else if (NUMCLUSTERS == 3)
                {
                    if (OUTPUTRT)
                    {
                        header1 = ("Scan number(s), Raw File(s), Charge State(s), Best Scan Number, Peptide, # Labels, Cluster 1 Avg Theo MZ, Cluster 2 Avg Theo MZ, Cluster 3 Avg Theo MZ, Cluster 1 Avg Adjusted Theo MZ, Cluster 2 Avg Adjusted Theo MZ, Cluster 3 Avg Adjusted Theo MZ, Cluster 1 Resolvable?, Cluster 2 Resolvable?, Cluster 3 Resolvable?, Quantified MS1 Scans, Cluster 1 # Measurements, Cluster 2 # Measurements, Cluster 3 # Measurements, Intensity 1, Intensity 2, Intensity 3, Intensity 4, Intensity 5, Intensity 6, RT 1, RT 2, RT 3, RT 4, RT 5, RT 6, Ratio 2/1, Ratio 4/3, Ratio 6/5, Cluster 1 Ratio Count, Cluster 2 Ratio Count, Cluster 3 Ratio Count, Cluster 1 Missing Channels?, Cluster 2 Missing Channels?, Cluster 3 Missing Channels?, Quant Reason");
                        writer1.WriteLine(header1);
                        foreach (PeptideID peptide in allPeptides.Values)
                        {
                            writer1.WriteLine("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23},{24},{25},{26},{27},{28},{29},{30},{31},{32},{33},{34},{35},{36},{37},{38},{39},{40}",
                                peptide.scanNumbers, peptide.rawFiles, peptide.chargeStates, peptide.bestPSM.ScanNumber, peptide.sequence, peptide.numLabels,
                                peptide.averageTheoMZ[0], peptide.averageTheoMZ[1], peptide.averageTheoMZ[2], peptide.averageAdjustedTheoMZ[0], peptide.averageAdjustedTheoMZ[1], peptide.averageAdjustedTheoMZ[2],
                                peptide.GetTheoreticalResolvability(0), peptide.GetTheoreticalResolvability(1), peptide.GetTheoreticalResolvability(2),
                                peptide.countAllPairs, peptide.countAllIsotopes[0], peptide.countAllIsotopes[1], peptide.countAllIsotopes[2], peptide.totalIntensity[0, NUMISOTOPES], peptide.totalIntensity[1, NUMISOTOPES], peptide.totalIntensity[2, NUMISOTOPES], peptide.totalIntensity[3, NUMISOTOPES], peptide.totalIntensity[4, NUMISOTOPES], peptide.totalIntensity[5, NUMISOTOPES],
                                peptide.maxCompleteIntensity[0, 1], peptide.maxCompleteIntensity[1, 1], peptide.maxCompleteIntensity[2, 1], peptide.maxCompleteIntensity[3, 1], peptide.maxCompleteIntensity[4, 1], peptide.maxCompleteIntensity[5, 1],
                                peptide.heavyToLightRatioSum[0, 0], peptide.heavyToLightRatioSum[1, 0], peptide.heavyToLightRatioSum[2, 0], peptide.finalQuantified[0], peptide.finalQuantified[1], peptide.finalQuantified[2], peptide.quantifiedNoiseIncluded[0], peptide.quantifiedNoiseIncluded[1], peptide.quantifiedNoiseIncluded[2], peptide.noQuantReason);
                        }
                    }
                    else
                    {
                        header1 = ("Scan number(s), Raw File(s), Charge State(s), Best Scan Number, Peptide, # Labels, Cluster 1 Avg Theo MZ, Cluster 2 Avg Theo MZ, Cluster 3 Avg Theo MZ, Cluster 1 Avg Adjusted Theo MZ, Cluster 2 Avg Adjusted Theo MZ, Cluster 3 Avg Adjusted Theo MZ, Cluster 1 Resolvable?, Cluster 2 Resolvable?, Cluster 3 Resolvable?, Quantified MS1 Scans, Cluster 1 # Measurements, Cluster 2 # Measurements, Cluster 3 # Measurements, Intensity 1, Intensity 2, Intensity 3, Intensity 4, Intensity 5, Intensity 6, Ratio 2/1, Ratio 4/3, Ratio 6/5, Cluster 1 Ratio Count, Cluster 2 Ratio Count, Cluster 3 Ratio Count, Cluster 1 Missing Channels?, Cluster 2 Missing Channels?, Cluster 3 Missing Channels?, Quant Reason");
                        writer1.WriteLine(header1);
                        foreach (PeptideID peptide in allPeptides.Values)
                        {
                            writer1.WriteLine("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23},{24},{25},{26},{27},{28},{29},{30},{31},{32},{33},{34}",
                                peptide.scanNumbers, peptide.rawFiles, peptide.chargeStates, peptide.bestPSM.ScanNumber, peptide.sequence, peptide.numLabels,
                                peptide.averageTheoMZ[0], peptide.averageTheoMZ[1], peptide.averageTheoMZ[2], peptide.averageAdjustedTheoMZ[0], peptide.averageAdjustedTheoMZ[1], peptide.averageAdjustedTheoMZ[2],
                                peptide.GetTheoreticalResolvability(0), peptide.GetTheoreticalResolvability(1), peptide.GetTheoreticalResolvability(2),
                                peptide.countAllPairs, peptide.countAllIsotopes[0], peptide.countAllIsotopes[1], peptide.countAllIsotopes[2], peptide.totalIntensity[0, NUMISOTOPES], peptide.totalIntensity[1, NUMISOTOPES], peptide.totalIntensity[2, NUMISOTOPES], peptide.totalIntensity[3, NUMISOTOPES], peptide.totalIntensity[4, NUMISOTOPES], peptide.totalIntensity[5, NUMISOTOPES],
                                peptide.heavyToLightRatioSum[0, 0], peptide.heavyToLightRatioSum[1, 0], peptide.heavyToLightRatioSum[2, 0], peptide.finalQuantified[0], peptide.finalQuantified[1], peptide.finalQuantified[2], peptide.quantifiedNoiseIncluded[0], peptide.quantifiedNoiseIncluded[1], peptide.quantifiedNoiseIncluded[2], peptide.noQuantReason);
                        }
                    }
                }
                writer1.Close();
            }
            else if (NUMCHANNELS == 8)
            {
                string header1;
                if (NUMCLUSTERS == 2)
                {
                    if (OUTPUTRT)
                    {
                        header1 = ("Scan number(s), Raw File(s), Charge State(s), Best Scan Number, Peptide, # Labels, Cluster 1 Avg Theo MZ, Cluster 2 Avg Theo MZ, Cluster 1 Avg Adjusted Theo MZ, Cluster 2 Avg Adjusted Theo MZ, Cluster 1 Resolvable?, Cluster 2 Resolvable?, Quantified MS1 Scans, Cluster 1 # Measurements, Cluster 2 # Measurements, Intensity 1, Intensity 2, Intensity 3, Intensity 4, Intensity 5, Intensity 6, Intensity 7, Intensity 8, RT 1, RT 2, RT 3, RT 4, RT 5, RT 6, RT 7, RT 8, Ratio 2/1, Ratio 3/1, Ratio 4/1, Ratio 6/5, Ratio 7/5, Ratio 8/5, Cluster 1 Ratio Count, Cluster 2 Ratio Count, Cluster 1 Missing Channels?, Cluster 2 Missing Channels?, Quant Reason");
                        writer1.WriteLine(header1);
                        foreach (PeptideID peptide in allPeptides.Values)
                        {
                            writer1.WriteLine("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23},{24},{25},{26},{27},{28},{29},{30},{31},{32},{33},{34},{35},{36},{37},{38},{39},{40},{41}",
                                peptide.scanNumbers, peptide.rawFiles, peptide.chargeStates, peptide.bestPSM.ScanNumber, peptide.sequence, peptide.numLabels,
                                peptide.averageTheoMZ[0], peptide.averageTheoMZ[1], peptide.averageAdjustedTheoMZ[0], peptide.averageAdjustedTheoMZ[1],
                                peptide.GetTheoreticalResolvability(0), peptide.GetTheoreticalResolvability(1),
                                peptide.countAllPairs, peptide.countAllIsotopes[0], peptide.countAllIsotopes[1], peptide.totalIntensity[0, NUMISOTOPES], peptide.totalIntensity[1, NUMISOTOPES], peptide.totalIntensity[2, NUMISOTOPES], peptide.totalIntensity[3, NUMISOTOPES], peptide.totalIntensity[4, NUMISOTOPES], peptide.totalIntensity[5, NUMISOTOPES], peptide.totalIntensity[6, NUMISOTOPES], peptide.totalIntensity[7, NUMISOTOPES],
                                peptide.maxCompleteIntensity[0, 1], peptide.maxCompleteIntensity[1, 1], peptide.maxCompleteIntensity[2, 1], peptide.maxCompleteIntensity[3, 1], peptide.maxCompleteIntensity[4, 1], peptide.maxCompleteIntensity[5, 1], peptide.maxCompleteIntensity[6,1], peptide.maxCompleteIntensity[7,1],
                                peptide.heavyToLightRatioSum[0, 0], peptide.heavyToLightRatioSum[1, 0], peptide.heavyToLightRatioSum[2, 0], peptide.heavyToLightRatioSum[3, 0], peptide.heavyToLightRatioSum[4,0], peptide.heavyToLightRatioSum[5,0], peptide.finalQuantified[0], peptide.finalQuantified[1], peptide.quantifiedNoiseIncluded[0], peptide.quantifiedNoiseIncluded[1], peptide.noQuantReason);
                        }
                    }
                    else
                    {
                        header1 = ("Scan number(s), Raw File(s), Charge State(s), Best Scan Number, Peptide, # Labels, Cluster 1 Avg Theo MZ, Cluster 2 Avg Theo MZ, Cluster 1 Avg Adjusted Theo MZ, Cluster 2 Avg Adjusted Theo MZ, Cluster 1 Resolvable?, Cluster 2 Resolvable?, Quantified MS1 Scans, Cluster 1 # Measurements, Cluster 2 # Measurements, Intensity 1, Intensity 2, Intensity 3, Intensity 4, Intensity 5, Intensity 6, Intensity 7, Intensity 8, Ratio 2/1, Ratio 3/1, Ratio 4/1, Ratio 6/5, Ratio 7/5, Ratio 8/5, Cluster 1 Ratio Count, Cluster 2 Ratio Count, Cluster 1 Missing Channels?, Cluster 2 Missing Channels?, Quant Reason");
                        writer1.WriteLine(header1);
                        foreach (PeptideID peptide in allPeptides.Values)
                        {
                            writer1.WriteLine("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23},{24},{25},{26},{27},{28},{29},{30},{31},{32},{33}",
                                peptide.scanNumbers, peptide.rawFiles, peptide.chargeStates, peptide.bestPSM.ScanNumber, peptide.sequence, peptide.numLabels,
                                peptide.averageTheoMZ[0], peptide.averageTheoMZ[1], peptide.averageAdjustedTheoMZ[0], peptide.averageAdjustedTheoMZ[1],
                                peptide.GetTheoreticalResolvability(0), peptide.GetTheoreticalResolvability(1),
                                peptide.countAllPairs, peptide.countAllIsotopes[0], peptide.countAllIsotopes[1], peptide.totalIntensity[0, NUMISOTOPES], peptide.totalIntensity[1, NUMISOTOPES], peptide.totalIntensity[2, NUMISOTOPES], peptide.totalIntensity[3, NUMISOTOPES], peptide.totalIntensity[4, NUMISOTOPES], peptide.totalIntensity[5, NUMISOTOPES], peptide.totalIntensity[6, NUMISOTOPES], peptide.totalIntensity[7, NUMISOTOPES],
                                peptide.heavyToLightRatioSum[0, 0], peptide.heavyToLightRatioSum[1, 0], peptide.heavyToLightRatioSum[2, 0], peptide.heavyToLightRatioSum[3, 0], peptide.heavyToLightRatioSum[4, 0], peptide.heavyToLightRatioSum[5, 0], peptide.finalQuantified[0], peptide.finalQuantified[1], peptide.quantifiedNoiseIncluded[0], peptide.quantifiedNoiseIncluded[1], peptide.noQuantReason);
                        }
                    }
                }
                writer1.Close();
            }
            else if (NUMCHANNELS == 9)
            {
                string header1;
                if (NUMCLUSTERS == 3)
                {
                    if (OUTPUTRT)
                    {
                        header1 = ("Scan number(s), Raw File(s), Charge State(s), Best Scan Number, Peptide, # Labels, Cluster 1 Avg Theo MZ, Cluster 2 Avg Theo MZ, Cluster 3 Avg Theo MZ, Cluster 1 Avg Adjusted Theo MZ, Cluster 2 Avg Adjusted Theo MZ, Cluster 3 Avg Adjusted Theo MZ, Cluster 1 Resolvable?, Cluster 2 Resolvable?, Cluster 3 Resolvable?, Quantified MS1 Scans, Cluster 1 # Measurements, Cluster 2 # Measurements, Cluster 3 # Measurements, Intensity 1, Intensity 2, Intensity 3, Intensity 4, Intensity 5, Intensity 6, Intensity 7, Intensity 8, Intensity 9, RT 1, RT 2, RT 3, RT 4, RT 5, RT 6, RT 7, RT 8, RT 9, Ratio 2/1, Ratio 3/1, Ratio 5/4, Ratio 6/4, Ratio 8/7, Ratio 9/7, Cluster 1 Ratio Count, Cluster 2 Ratio Count, Cluster 3 Ratio Count, Cluster 1 Missing Channels?, Cluster 2 Missing Channels?, Cluster 3 Missing Channels?, Quant Reason");
                        writer1.WriteLine(header1);
                        foreach (PeptideID peptide in allPeptides.Values)
                        {
                            writer1.WriteLine("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23},{24},{25},{26},{27},{28},{29},{30},{31},{32},{33},{34},{35},{36},{37},{38},{39},{40},{41},{42},{43},{44},{45},{46},{47},{48},{49}",
                                peptide.scanNumbers, peptide.rawFiles, peptide.chargeStates, peptide.bestPSM.ScanNumber, peptide.sequence, peptide.numLabels,
                                peptide.averageTheoMZ[0], peptide.averageTheoMZ[1], peptide.averageTheoMZ[2], peptide.averageAdjustedTheoMZ[0], peptide.averageAdjustedTheoMZ[1], peptide.averageAdjustedTheoMZ[2],
                                peptide.GetTheoreticalResolvability(0), peptide.GetTheoreticalResolvability(1), peptide.GetTheoreticalResolvability(2),
                                peptide.countAllPairs, peptide.countAllIsotopes[0], peptide.countAllIsotopes[1], peptide.countAllIsotopes[2], peptide.totalIntensity[0, NUMISOTOPES], peptide.totalIntensity[1, NUMISOTOPES], peptide.totalIntensity[2, NUMISOTOPES], peptide.totalIntensity[3, NUMISOTOPES], peptide.totalIntensity[4, NUMISOTOPES], peptide.totalIntensity[5, NUMISOTOPES], peptide.totalIntensity[6, NUMISOTOPES], peptide.totalIntensity[7, NUMISOTOPES], peptide.totalIntensity[8, NUMISOTOPES],
                                peptide.maxCompleteIntensity[0, 1], peptide.maxCompleteIntensity[1, 1], peptide.maxCompleteIntensity[2, 1], peptide.maxCompleteIntensity[3, 1], peptide.maxCompleteIntensity[4, 1], peptide.maxCompleteIntensity[5, 1], peptide.maxCompleteIntensity[6, 1], peptide.maxCompleteIntensity[7, 1], peptide.maxCompleteIntensity[8, 1],
                                peptide.heavyToLightRatioSum[0, 0], peptide.heavyToLightRatioSum[1, 0], peptide.heavyToLightRatioSum[2, 0], peptide.heavyToLightRatioSum[3,0], peptide.heavyToLightRatioSum[4,0], peptide.heavyToLightRatioSum[5,0], peptide.finalQuantified[0], peptide.finalQuantified[1], peptide.finalQuantified[2], peptide.quantifiedNoiseIncluded[0], peptide.quantifiedNoiseIncluded[1], peptide.quantifiedNoiseIncluded[2], peptide.noQuantReason);
                        }
                    }
                    else
                    {
                        header1 = ("Scan number(s), Raw File(s), Charge State(s), Best Scan Number, Peptide, # Labels, Cluster 1 Avg Theo MZ, Cluster 2 Avg Theo MZ, Cluster 3 Avg Theo MZ, Cluster 1 Avg Adjusted Theo MZ, Cluster 2 Avg Adjusted Theo MZ, Cluster 3 Avg Adjusted Theo MZ, Cluster 1 Resolvable?, Cluster 2 Resolvable?, Cluster 3 Resolvable?, Quantified MS1 Scans, Cluster 1 # Measurements, Cluster 2 # Measurements, Cluster 3 # Measurements, Intensity 1, Intensity 2, Intensity 3, Intensity 4, Intensity 5, Intensity 6, Intensity 7, Intensity 8, Intensity 9, Ratio 2/1, Ratio 3/1, Ratio 5/4, Ratio 6/4, Ratio 8/7, Ratio 9/7, Cluster 1 Ratio Count, Cluster 2 Ratio Count, Cluster 3 Ratio Count, Cluster 1 Missing Channels?, Cluster 2 Missing Channels?, Cluster 3 Missing Channels?, Quant Reason");
                        writer1.WriteLine(header1);
                        foreach (PeptideID peptide in allPeptides.Values)
                        {
                            writer1.WriteLine("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23},{24},{25},{26},{27},{28},{29},{30},{31},{32},{33},{34},{35},{36},{37},{38},{39},{40}",
                                peptide.scanNumbers, peptide.rawFiles, peptide.chargeStates, peptide.bestPSM.ScanNumber, peptide.sequence, peptide.numLabels,
                                peptide.averageTheoMZ[0], peptide.averageTheoMZ[1], peptide.averageTheoMZ[2], peptide.averageAdjustedTheoMZ[0], peptide.averageAdjustedTheoMZ[1], peptide.averageAdjustedTheoMZ[2],
                                peptide.GetTheoreticalResolvability(0), peptide.GetTheoreticalResolvability(1), peptide.GetTheoreticalResolvability(2),
                                peptide.countAllPairs, peptide.countAllIsotopes[0], peptide.countAllIsotopes[1], peptide.countAllIsotopes[2], peptide.totalIntensity[0, NUMISOTOPES], peptide.totalIntensity[1, NUMISOTOPES], peptide.totalIntensity[2, NUMISOTOPES], peptide.totalIntensity[3, NUMISOTOPES], peptide.totalIntensity[4, NUMISOTOPES], peptide.totalIntensity[5, NUMISOTOPES], peptide.totalIntensity[6, NUMISOTOPES], peptide.totalIntensity[7, NUMISOTOPES], peptide.totalIntensity[8, NUMISOTOPES],
                                peptide.heavyToLightRatioSum[0, 0], peptide.heavyToLightRatioSum[1, 0], peptide.heavyToLightRatioSum[2, 0], peptide.heavyToLightRatioSum[3, 0], peptide.heavyToLightRatioSum[4, 0], peptide.heavyToLightRatioSum[5, 0], peptide.finalQuantified[0], peptide.finalQuantified[1], peptide.finalQuantified[2], peptide.quantifiedNoiseIncluded[0], peptide.quantifiedNoiseIncluded[1], peptide.quantifiedNoiseIncluded[2], peptide.noQuantReason);
                        }
                    }
                }
                writer1.Close();
            }
            else if (NUMCHANNELS == 12)
            {
                string header1;
                if (NUMCLUSTERS == 2)
                {
                    if (OUTPUTRT)
                    {
                        header1 = ("Scan number(s), Raw File(s), Charge State(s), Best Scan Number, Peptide, # Labels, Cluster 1 Avg Theo MZ, Cluster 2 Avg Theo MZ, Cluster 1 Avg Adjusted Theo MZ, Cluster 2 Avg Adjusted Theo MZ, Cluster 1 Resolvable?, Cluster 2 Resolvable?, Quantified MS1 Scans, Cluster 1 # Measurements, Cluster 2 # Measurements, Intensity 1, Intensity 2, Intensity 3, Intensity 4, Intensity 5, Intensity 6, Intensity 7, Intensity 8, Intensity 9, Intensity 10, Intensity 11, Intensity 12, RT 1, RT 2, RT 3, RT 4, RT 5, RT 6, RT 7, RT 8, RT 9, RT 10, RT 11, RT 12, Ratio 2/1, Ratio 3/1, Ratio 4/1, Ratio 5/1, Ratio 6/1, Ratio 8/7, Ratio 9/7, Ratio 10/7, Ratio 11/7, Ratio 12/7, Cluster 1 Ratio Count, Cluster 2 Ratio Count, Cluster 1 Missing Channels?, Cluster 2 Missing Channels?, Quant Reason");
                        writer1.WriteLine(header1);
                        foreach (PeptideID peptide in allPeptides.Values)
                        {
                            writer1.WriteLine("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23},{24},{25},{26},{27},{28},{29},{30},{31},{32},{33},{34},{35},{36},{37},{38},{39},{40},{41},{42},{43},{44},{45},{46},{47},{48},{49},{50},{51},{52},{53}",
                                peptide.scanNumbers, peptide.rawFiles, peptide.chargeStates, peptide.bestPSM.ScanNumber, peptide.sequence, peptide.numLabels,
                                peptide.averageTheoMZ[0], peptide.averageTheoMZ[1], peptide.averageAdjustedTheoMZ[0], peptide.averageAdjustedTheoMZ[1],
                                peptide.GetTheoreticalResolvability(0), peptide.GetTheoreticalResolvability(1),
                                peptide.countAllPairs, peptide.countAllIsotopes[0], peptide.countAllIsotopes[1], peptide.totalIntensity[0, NUMISOTOPES], peptide.totalIntensity[1, NUMISOTOPES], peptide.totalIntensity[2, NUMISOTOPES], peptide.totalIntensity[3, NUMISOTOPES], peptide.totalIntensity[4, NUMISOTOPES], peptide.totalIntensity[5, NUMISOTOPES], peptide.totalIntensity[6, NUMISOTOPES], peptide.totalIntensity[7, NUMISOTOPES], peptide.totalIntensity[8, NUMISOTOPES], peptide.totalIntensity[9, NUMISOTOPES], peptide.totalIntensity[10, NUMISOTOPES], peptide.totalIntensity[11, NUMISOTOPES],
                                peptide.maxCompleteIntensity[0, 1], peptide.maxCompleteIntensity[1, 1], peptide.maxCompleteIntensity[2, 1], peptide.maxCompleteIntensity[3, 1], peptide.maxCompleteIntensity[4, 1], peptide.maxCompleteIntensity[5, 1], peptide.maxCompleteIntensity[6, 1], peptide.maxCompleteIntensity[7, 1], peptide.maxCompleteIntensity[8,1], peptide.maxCompleteIntensity[9,1], peptide.maxCompleteIntensity[10,1], peptide.maxCompleteIntensity[11,1],
                                peptide.heavyToLightRatioSum[0, 0], peptide.heavyToLightRatioSum[1, 0], peptide.heavyToLightRatioSum[2, 0], peptide.heavyToLightRatioSum[3, 0], peptide.heavyToLightRatioSum[4, 0], peptide.heavyToLightRatioSum[5, 0], peptide.heavyToLightRatioSum[6,0], peptide.heavyToLightRatioSum[7,0], peptide.heavyToLightRatioSum[8,0], peptide.heavyToLightRatioSum[9,0], peptide.finalQuantified[0], peptide.finalQuantified[1], peptide.quantifiedNoiseIncluded[0], peptide.quantifiedNoiseIncluded[1], peptide.noQuantReason);
                        }
                    }
                    else
                    {
                        header1 = ("Scan number(s), Raw File(s), Charge State(s), Best Scan Number, Peptide, # Labels, Cluster 1 Avg Theo MZ, Cluster 2 Avg Theo MZ, Cluster 1 Avg Adjusted Theo MZ, Cluster 2 Avg Adjusted Theo MZ, Cluster 1 Resolvable?, Cluster 2 Resolvable?, Quantified MS1 Scans, Cluster 1 # Measurements, Cluster 2 # Measurements, Intensity 1, Intensity 2, Intensity 3, Intensity 4, Intensity 5, Intensity 6, Intensity 7, Intensity 8, Intensity 9, Intensity 10, Intensity 11, Intensity 12, Ratio 2/1, Ratio 3/1, Ratio 4/1, Ratio 5/1, Ratio 6/1, Ratio 8/7, Ratio 9/7, Ratio 10/7, Ratio 11/7, Ratio 12/7, Cluster 1 Ratio Count, Cluster 2 Ratio Count, Cluster 1 Missing Channels?, Cluster 2 Missing Channels?, Quant Reason");
                        writer1.WriteLine(header1);
                        foreach (PeptideID peptide in allPeptides.Values)
                        {
                            writer1.WriteLine("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23},{24},{25},{26},{27},{28},{29},{30},{31},{32},{33},{34},{35},{36},{37},{38},{39},{40},{41}",
                                peptide.scanNumbers, peptide.rawFiles, peptide.chargeStates, peptide.bestPSM.ScanNumber, peptide.sequence, peptide.numLabels,
                                peptide.averageTheoMZ[0], peptide.averageTheoMZ[1], peptide.averageAdjustedTheoMZ[0], peptide.averageAdjustedTheoMZ[1],
                                peptide.GetTheoreticalResolvability(0), peptide.GetTheoreticalResolvability(1),
                                peptide.countAllPairs, peptide.countAllIsotopes[0], peptide.countAllIsotopes[1], peptide.totalIntensity[0, NUMISOTOPES], peptide.totalIntensity[1, NUMISOTOPES], peptide.totalIntensity[2, NUMISOTOPES], peptide.totalIntensity[3, NUMISOTOPES], peptide.totalIntensity[4, NUMISOTOPES], peptide.totalIntensity[5, NUMISOTOPES], peptide.totalIntensity[6, NUMISOTOPES], peptide.totalIntensity[7, NUMISOTOPES], peptide.totalIntensity[8, NUMISOTOPES], peptide.totalIntensity[9, NUMISOTOPES], peptide.totalIntensity[10, NUMISOTOPES], peptide.totalIntensity[11, NUMISOTOPES],
                                peptide.heavyToLightRatioSum[0, 0], peptide.heavyToLightRatioSum[1, 0], peptide.heavyToLightRatioSum[2, 0], peptide.heavyToLightRatioSum[3, 0], peptide.heavyToLightRatioSum[4, 0], peptide.heavyToLightRatioSum[5, 0], peptide.heavyToLightRatioSum[6, 0], peptide.heavyToLightRatioSum[7, 0], peptide.heavyToLightRatioSum[8, 0], peptide.heavyToLightRatioSum[9, 0], peptide.finalQuantified[0], peptide.finalQuantified[1], peptide.quantifiedNoiseIncluded[0], peptide.quantifiedNoiseIncluded[1], peptide.noQuantReason);
                        }
                    }
                }
                else if (NUMCLUSTERS == 3)
                {
                    if (OUTPUTRT)
                    {
                        header1 = ("Scan number(s), Raw File(s), Charge State(s), Best Scan Number, Peptide, # Labels, Cluster 1 Avg Theo MZ, Cluster 2 Avg Theo MZ, Cluster 3 Avg Theo MZ, Cluster 1 Avg Adjusted Theo MZ, Cluster 2 Avg Adjusted Theo MZ, Cluster 3 Avg Adjusted Theo MZ, Cluster 1 Resolvable?, Cluster 2 Resolvable?, Cluster 3 Resolvable?, Quantified MS1 Scans, Cluster 1 # Measurements, Cluster 2 # Measurements, Cluster 3 # Measurements, Intensity 1, Intensity 2, Intensity 3, Intensity 4, Intensity 5, Intensity 6, Intensity 7, Intensity 8, Intensity 9, Intensity 10, Intensity 11, Intensity 12, RT 1, RT 2, RT 3, RT 4, RT 5, RT 6, RT 7, RT 8, RT 9, RT 10, RT 11, RT 12, Ratio 2/1, Ratio 3/1, Ratio 4/1, Ratio 6/5, Ratio 7/5, Ratio 8/5, Ratio 10/9, Ratio 11/9, Ratio 12/9, Cluster 1 Ratio Count, Cluster 2 Ratio Count, Cluster 3 Ratio Count, Cluster 1 Missing Channels?, Cluster 2 Missing Channels?, Cluster 3 Missing Channels?, Quant Reason");
                        writer1.WriteLine(header1);
                        foreach (PeptideID peptide in allPeptides.Values)
                        {
                            writer1.WriteLine("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23},{24},{25},{26},{27},{28},{29},{30},{31},{32},{33},{34},{35},{36},{37},{38},{39},{40},{41},{42},{43},{44},{45},{46},{47},{48},{49},{50},{51},{52},{53},{54},{55},{56},{57},{58}",
                                peptide.scanNumbers, peptide.rawFiles, peptide.chargeStates, peptide.bestPSM.ScanNumber, peptide.sequence, peptide.numLabels,
                                peptide.averageTheoMZ[0], peptide.averageTheoMZ[1], peptide.averageTheoMZ[2], peptide.averageAdjustedTheoMZ[0], peptide.averageAdjustedTheoMZ[1], peptide.averageAdjustedTheoMZ[2],
                                peptide.GetTheoreticalResolvability(0), peptide.GetTheoreticalResolvability(1), peptide.GetTheoreticalResolvability(2),
                                peptide.countAllPairs, peptide.countAllIsotopes[0], peptide.countAllIsotopes[1], peptide.countAllIsotopes[2], peptide.totalIntensity[0, NUMISOTOPES], peptide.totalIntensity[1, NUMISOTOPES], peptide.totalIntensity[2, NUMISOTOPES], peptide.totalIntensity[3, NUMISOTOPES], peptide.totalIntensity[4, NUMISOTOPES], peptide.totalIntensity[5, NUMISOTOPES], peptide.totalIntensity[6, NUMISOTOPES], peptide.totalIntensity[7, NUMISOTOPES], peptide.totalIntensity[8, NUMISOTOPES], peptide.totalIntensity[9, NUMISOTOPES], peptide.totalIntensity[10, NUMISOTOPES], peptide.totalIntensity[11, NUMISOTOPES],
                                peptide.maxCompleteIntensity[0, 1], peptide.maxCompleteIntensity[1, 1], peptide.maxCompleteIntensity[2, 1], peptide.maxCompleteIntensity[3, 1], peptide.maxCompleteIntensity[4, 1], peptide.maxCompleteIntensity[5, 1], peptide.maxCompleteIntensity[6, 1], peptide.maxCompleteIntensity[7, 1], peptide.maxCompleteIntensity[8, 1], peptide.maxCompleteIntensity[9,1], peptide.maxCompleteIntensity[10,1], peptide.maxCompleteIntensity[11,1],
                                peptide.heavyToLightRatioSum[0, 0], peptide.heavyToLightRatioSum[1, 0], peptide.heavyToLightRatioSum[2, 0], peptide.heavyToLightRatioSum[3, 0], peptide.heavyToLightRatioSum[4, 0], peptide.heavyToLightRatioSum[5, 0], peptide.heavyToLightRatioSum[6,0], peptide.heavyToLightRatioSum[7,0], peptide.heavyToLightRatioSum[8,0], peptide.finalQuantified[0], peptide.finalQuantified[1], peptide.finalQuantified[2], peptide.quantifiedNoiseIncluded[0], peptide.quantifiedNoiseIncluded[1], peptide.quantifiedNoiseIncluded[2], peptide.noQuantReason);
                        }
                    }
                    else
                    {
                        header1 = ("Scan number(s), Raw File(s), Charge State(s), Best Scan Number, Peptide, # Labels, Cluster 1 Avg Theo MZ, Cluster 2 Avg Theo MZ, Cluster 3 Avg Theo MZ, Cluster 1 Avg Adjusted Theo MZ, Cluster 2 Avg Adjusted Theo MZ, Cluster 3 Avg Adjusted Theo MZ, Cluster 1 Resolvable?, Cluster 2 Resolvable?, Cluster 3 Resolvable?, Quantified MS1 Scans, Cluster 1 # Measurements, Cluster 2 # Measurements, Cluster 3 # Measurements, Intensity 1, Intensity 2, Intensity 3, Intensity 4, Intensity 5, Intensity 6, Intensity 7, Intensity 8, Intensity 9, Intensity 10, Intensity 11, Intensity 12, Ratio 2/1, Ratio 3/1, Ratio 4/1, Ratio 6/5, Ratio 7/5, Ratio 8/5, Ratio 10/9, Ratio 11/9, Ratio 12/9, Cluster 1 Ratio Count, Cluster 2 Ratio Count, Cluster 3 Ratio Count, Cluster 1 Missing Channels?, Cluster 2 Missing Channels?, Cluster 3 Missing Channels?, Quant Reason");
                        writer1.WriteLine(header1);
                        foreach (PeptideID peptide in allPeptides.Values)
                        {
                            writer1.WriteLine("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23},{24},{25},{26},{27},{28},{29},{30},{31},{32},{33},{34},{35},{36},{37},{38},{39},{40},{41},{42},{43},{44},{45},{46}",
                                peptide.scanNumbers, peptide.rawFiles, peptide.chargeStates, peptide.bestPSM.ScanNumber, peptide.sequence, peptide.numLabels,
                                peptide.averageTheoMZ[0], peptide.averageTheoMZ[1], peptide.averageTheoMZ[2], peptide.averageAdjustedTheoMZ[0], peptide.averageAdjustedTheoMZ[1], peptide.averageAdjustedTheoMZ[2],
                                peptide.GetTheoreticalResolvability(0), peptide.GetTheoreticalResolvability(1), peptide.GetTheoreticalResolvability(2),
                                peptide.countAllPairs, peptide.countAllIsotopes[0], peptide.countAllIsotopes[1], peptide.countAllIsotopes[2], peptide.totalIntensity[0, NUMISOTOPES], peptide.totalIntensity[1, NUMISOTOPES], peptide.totalIntensity[2, NUMISOTOPES], peptide.totalIntensity[3, NUMISOTOPES], peptide.totalIntensity[4, NUMISOTOPES], peptide.totalIntensity[5, NUMISOTOPES], peptide.totalIntensity[6, NUMISOTOPES], peptide.totalIntensity[7, NUMISOTOPES], peptide.totalIntensity[8, NUMISOTOPES], peptide.totalIntensity[9, NUMISOTOPES], peptide.totalIntensity[10, NUMISOTOPES], peptide.totalIntensity[11, NUMISOTOPES],
                                peptide.heavyToLightRatioSum[0, 0], peptide.heavyToLightRatioSum[1, 0], peptide.heavyToLightRatioSum[2, 0], peptide.heavyToLightRatioSum[3, 0], peptide.heavyToLightRatioSum[4, 0], peptide.heavyToLightRatioSum[5, 0], peptide.heavyToLightRatioSum[6, 0], peptide.heavyToLightRatioSum[7, 0], peptide.heavyToLightRatioSum[8, 0], peptide.finalQuantified[0], peptide.finalQuantified[1], peptide.finalQuantified[2], peptide.quantifiedNoiseIncluded[0], peptide.quantifiedNoiseIncluded[1], peptide.quantifiedNoiseIncluded[2], peptide.noQuantReason);
                        }
                    }
                }
                writer1.Close();
            }
            else if (NUMCHANNELS == 18)
            {
                string header1;
                if (NUMCLUSTERS == 3)
                {
                    if (OUTPUTRT)
                    {
                        header1 = ("Scan number(s), Raw File(s), Charge State(s), Best Scan Number, Peptide, # Labels, Cluster 1 Avg Theo MZ, Cluster 2 Avg Theo MZ, Cluster 3 Avg Theo MZ, Cluster 1 Avg Adjusted Theo MZ, Cluster 2 Avg Adjusted Theo MZ, Cluster 3 Avg Adjusted Theo MZ, Cluster 1 Resolvable?, Cluster 2 Resolvable?, Cluster 3 Resolvable?, Quantified MS1 Scans, Cluster 1 # Measurements, Cluster 2 # Measurements, Cluster 3 # Measurements, Intensity 1, Intensity 2, Intensity 3, Intensity 4, Intensity 5, Intensity 6, Intensity 7, Intensity 8, Intensity 9, Intensity 10, Intensity 11, Intensity 12, Intensity 13, Intensity 14, Intensity 15, Intensity 16, Intensity 17, Intensity 18, RT 1, RT 2, RT 3, RT 4, RT 5, RT 6, RT 7, RT 8, RT 9, RT 10, RT 11, RT 12, RT 13, RT 14, RT 15, RT 16, RT 17, RT 18, Ratio 2/1, Ratio 3/1, Ratio 4/1, Ratio 5/1, Ratio 6/1, Ratio 8/7, Ratio 9/7, Ratio 10/7, Ratio 11/7, Ratio 12/7, Ratio 14/13, Ratio 15/13, Ratio 16/13, Ratio 17/13, Ratio 18/13, Cluster 1 Ratio Count, Cluster 2 Ratio Count, Cluster 3 Ratio Count, Cluster 1 Missing Channels?, Cluster 2 Missing Channels?, Cluster 3 Missing Channels?, Quant Reason");
                        writer1.WriteLine(header1);
                        foreach (PeptideID peptide in allPeptides.Values)
                        {
                            writer1.WriteLine("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23},{24},{25},{26},{27},{28},{29},{30},{31},{32},{33},{34},{35},{36},{37},{38},{39},{40},{41},{42},{43},{44},{45},{46},{47},{48},{49},{50},{51},{52},{53},{54},{55},{56},{57},{58},{59},{60},{61},{62},{63},{64},{65},{66},{67,{68},{69},{70},{71},{72},{73},{74},{75},{76}",
                                peptide.scanNumbers, peptide.rawFiles, peptide.chargeStates, peptide.bestPSM.ScanNumber, peptide.sequence, peptide.numLabels,
                                peptide.averageTheoMZ[0], peptide.averageTheoMZ[1], peptide.averageTheoMZ[2], peptide.averageAdjustedTheoMZ[0], peptide.averageAdjustedTheoMZ[1], peptide.averageAdjustedTheoMZ[2],
                                peptide.GetTheoreticalResolvability(0), peptide.GetTheoreticalResolvability(1), peptide.GetTheoreticalResolvability(2),
                                peptide.countAllPairs, peptide.countAllIsotopes[0], peptide.countAllIsotopes[1], peptide.countAllIsotopes[2], peptide.totalIntensity[0, NUMISOTOPES], peptide.totalIntensity[1, NUMISOTOPES], peptide.totalIntensity[2, NUMISOTOPES], peptide.totalIntensity[3, NUMISOTOPES], peptide.totalIntensity[4, NUMISOTOPES], peptide.totalIntensity[5, NUMISOTOPES], peptide.totalIntensity[6, NUMISOTOPES], peptide.totalIntensity[7, NUMISOTOPES], peptide.totalIntensity[8, NUMISOTOPES], peptide.totalIntensity[9, NUMISOTOPES], peptide.totalIntensity[10, NUMISOTOPES], peptide.totalIntensity[11, NUMISOTOPES], peptide.totalIntensity[12, NUMISOTOPES], peptide.totalIntensity[13, NUMISOTOPES], peptide.totalIntensity[14, NUMISOTOPES], peptide.totalIntensity[15, NUMISOTOPES], peptide.totalIntensity[16, NUMISOTOPES], peptide.totalIntensity[17, NUMISOTOPES],
                                peptide.maxCompleteIntensity[0, 1], peptide.maxCompleteIntensity[1, 1], peptide.maxCompleteIntensity[2, 1], peptide.maxCompleteIntensity[3, 1], peptide.maxCompleteIntensity[4, 1], peptide.maxCompleteIntensity[5, 1], peptide.maxCompleteIntensity[6, 1], peptide.maxCompleteIntensity[7, 1], peptide.maxCompleteIntensity[8, 1], peptide.maxCompleteIntensity[9, 1], peptide.maxCompleteIntensity[10, 1], peptide.maxCompleteIntensity[11, 1], peptide.maxCompleteIntensity[12,1], peptide.maxCompleteIntensity[13,1], peptide.maxCompleteIntensity[14,1], peptide.maxCompleteIntensity[15,1], peptide.maxCompleteIntensity[16,1], peptide.maxCompleteIntensity[17,1],
                                peptide.heavyToLightRatioSum[0, 0], peptide.heavyToLightRatioSum[1, 0], peptide.heavyToLightRatioSum[2, 0], peptide.heavyToLightRatioSum[3, 0], peptide.heavyToLightRatioSum[4, 0], peptide.heavyToLightRatioSum[5, 0], peptide.heavyToLightRatioSum[6, 0], peptide.heavyToLightRatioSum[7, 0], peptide.heavyToLightRatioSum[8, 0], peptide.heavyToLightRatioSum[9,0], peptide.heavyToLightRatioSum[10,0], peptide.heavyToLightRatioSum[11,0], peptide.heavyToLightRatioSum[12,0], peptide.heavyToLightRatioSum[13,0], peptide.heavyToLightRatioSum[14,0], peptide.finalQuantified[0], peptide.finalQuantified[1], peptide.finalQuantified[2], peptide.quantifiedNoiseIncluded[0], peptide.quantifiedNoiseIncluded[1], peptide.quantifiedNoiseIncluded[2], peptide.noQuantReason);
                        }
                    }
                    else
                    {
                        header1 = ("Scan number(s), Raw File(s), Charge State(s), Best Scan Number, Peptide, # Labels, Cluster 1 Avg Theo MZ, Cluster 2 Avg Theo MZ, Cluster 3 Avg Theo MZ, Cluster 1 Avg Adjusted Theo MZ, Cluster 2 Avg Adjusted Theo MZ, Cluster 3 Avg Adjusted Theo MZ, Cluster 1 Resolvable?, Cluster 2 Resolvable?, Cluster 3 Resolvable?, Quantified MS1 Scans, Cluster 1 # Measurements, Cluster 2 # Measurements, Cluster 3 # Measurements, Intensity 1, Intensity 2, Intensity 3, Intensity 4, Intensity 5, Intensity 6, Intensity 7, Intensity 8, Intensity 9, Intensity 10, Intensity 11, Intensity 12, Intensity 13, Intensity 14, Intensity 15, Intensity 16, Intensity 17, Intensity 18, Ratio 2/1, Ratio 3/1, Ratio 4/1, Ratio 5/1, Ratio 6/1, Ratio 8/7, Ratio 9/7, Ratio 10/7, Ratio 11/7, Ratio 12/7, Ratio 14/13, Ratio 15/13, Ratio 16/13, Ratio 17/13, Ratio 18/13, Cluster 1 Ratio Count, Cluster 2 Ratio Count, Cluster 3 Ratio Count, Cluster 1 Missing Channels?, Cluster 2 Missing Channels?, Cluster 3 Missing Channels?, Quant Reason");
                        writer1.WriteLine(header1);
                        foreach (PeptideID peptide in allPeptides.Values)
                        {
                            writer1.WriteLine("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23},{24},{25},{26},{27},{28},{29},{30},{31},{32},{33},{34},{35},{36},{37},{38},{39},{40},{41},{42},{43},{44},{45},{46},{47},{48},{49},{50},{51},{52},{53},{54},{55},{56},{57},{58}",
                                peptide.scanNumbers, peptide.rawFiles, peptide.chargeStates, peptide.bestPSM.ScanNumber, peptide.sequence, peptide.numLabels,
                                peptide.averageTheoMZ[0], peptide.averageTheoMZ[1], peptide.averageTheoMZ[2], peptide.averageAdjustedTheoMZ[0], peptide.averageAdjustedTheoMZ[1], peptide.averageAdjustedTheoMZ[2],
                                peptide.GetTheoreticalResolvability(0), peptide.GetTheoreticalResolvability(1), peptide.GetTheoreticalResolvability(2),
                                peptide.countAllPairs, peptide.countAllIsotopes[0], peptide.countAllIsotopes[1], peptide.countAllIsotopes[2], peptide.totalIntensity[0, NUMISOTOPES], peptide.totalIntensity[1, NUMISOTOPES], peptide.totalIntensity[2, NUMISOTOPES], peptide.totalIntensity[3, NUMISOTOPES], peptide.totalIntensity[4, NUMISOTOPES], peptide.totalIntensity[5, NUMISOTOPES], peptide.totalIntensity[6, NUMISOTOPES], peptide.totalIntensity[7, NUMISOTOPES], peptide.totalIntensity[8, NUMISOTOPES], peptide.totalIntensity[9, NUMISOTOPES], peptide.totalIntensity[10, NUMISOTOPES], peptide.totalIntensity[11, NUMISOTOPES], peptide.totalIntensity[12, NUMISOTOPES], peptide.totalIntensity[13, NUMISOTOPES], peptide.totalIntensity[14, NUMISOTOPES], peptide.totalIntensity[15, NUMISOTOPES], peptide.totalIntensity[16, NUMISOTOPES], peptide.totalIntensity[17, NUMISOTOPES],
                                peptide.heavyToLightRatioSum[0, 0], peptide.heavyToLightRatioSum[1, 0], peptide.heavyToLightRatioSum[2, 0], peptide.heavyToLightRatioSum[3, 0], peptide.heavyToLightRatioSum[4, 0], peptide.heavyToLightRatioSum[5, 0], peptide.heavyToLightRatioSum[6, 0], peptide.heavyToLightRatioSum[7, 0], peptide.heavyToLightRatioSum[8, 0], peptide.heavyToLightRatioSum[9, 0], peptide.heavyToLightRatioSum[10, 0], peptide.heavyToLightRatioSum[11, 0], peptide.heavyToLightRatioSum[12, 0], peptide.heavyToLightRatioSum[13, 0], peptide.heavyToLightRatioSum[14, 0], peptide.finalQuantified[0], peptide.finalQuantified[1], peptide.finalQuantified[2], peptide.quantifiedNoiseIncluded[0], peptide.quantifiedNoiseIncluded[1], peptide.quantifiedNoiseIncluded[2], peptide.noQuantReason);
                        }
                    }
                }
                writer1.Close();
            }
        }

        private void writeTagQuantOutputFile(Dictionary<string, PeptideID> allPeptides, string file)
        {
            CsvReader psmFileReader = new CsvReader(new StreamReader(file), true);
            string tagQuantOutputFile = Path.Combine(outputFolderBox.Text, Path.GetFileNameWithoutExtension(file) + "_TagQuant.csv");
            StreamWriter tagQuantWriter = new StreamWriter(tagQuantOutputFile);

            string tagQuantHeader = ("Spectrum Number, Filename/id, Peptide, E-value, Mass, gi, Accession, Start, Stop, Defline, Mods, Charge, Theo Mass, P-value, NIST score, Precursor Isolation m/z, Precursor Isolation Mass (Da), Precursor Theoretical Neutral Mass (Da), Precursor Experimental Neutral Mass (Da), Precursor Mass Error (ppm), Adjusted Precursor Mass Error (ppm), Q-Value (%), Interference, TQ_126_1_NL, TQ_127_2_NL, TQ_128_3_NL, TQ_129_4_NL, TQ_130_5_NL, TQ_131_6_NL, TQ_126_1_dNL, TQ_127_2_dNL, TQ_128_3_dNL, TQ_129_4_dNL, TQ_130_5_dNL, TQ_131_6_dNL, TQ_Total_int(all_tags), TQ_126_1_PC, TQ_127_2_PC, TQ_128_3_PC, TQ_129_4_PC, TQ_130_5_PC, TQ_131_6_PC, TQ_126_1_PCN, TQ_127_2_PCN, TQ_128_3_PCN, TQ_129_4_PCN, TQ_130_5_PCN, TQ_131_6_PCN, Cleavage Site");
            tagQuantWriter.WriteLine(tagQuantHeader);

            bool includeInTQ;
            int charge;
            int scanNumber;
            string sequence;
            PeptideID quantifiedPeptide;

            while (psmFileReader.ReadNextRecord())
            {
                quantifiedPeptide = null;
                includeInTQ = false;
                sequence = psmFileReader["Peptide"].ToUpper();
                charge = int.Parse(psmFileReader["Charge"]);
                scanNumber = int.Parse(psmFileReader["Spectrum Number"]);

                foreach (PeptideID peptide in allPeptides.Values)
                {
                    if (!includeInTQ && peptide.sequence.ToUpper().Equals(sequence))
                    {
                        if (peptide.bestPSM.Charge == charge && peptide.bestPSM.ScanNumber == scanNumber)
                        {
                            includeInTQ = true;
                            quantifiedPeptide = peptide;
                        }
                    }
                }

                if (quantifiedPeptide != null)
                {
                    string defLineWithCommas = psmFileReader["Defline"];
                    if (defLineWithCommas.Contains(','))
                    {
                        defLineWithCommas = psmFileReader["Defline"].Replace(',', ';');
                    }

                    string modsWithCommas = psmFileReader["Mods"];
                    if (modsWithCommas.Contains(','))
                    {
                        modsWithCommas = psmFileReader["Mods"].Replace(',', ';');
                    }

                    if (NUMCLUSTERS == 1)
                    {
                        switch (NUMISOTOPOLOGUES)
                        {
                            case 2:
                                tagQuantWriter.WriteLine("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23},{24},{25},{26},{27},{28},{29},{30},{31},{32},{33},{34},{35},{36},{37},{38},{39},{40},{41},{42},{43},{44},{45},{46},{47},{48}",
                                    psmFileReader["Spectrum Number"], psmFileReader["Filename/id"], psmFileReader["Peptide"], psmFileReader["E-value"], psmFileReader["Mass"], psmFileReader["gi"], psmFileReader["Accession"], psmFileReader["Start"], psmFileReader["Stop"], defLineWithCommas, modsWithCommas,
                                    psmFileReader["Charge"], psmFileReader["Theo Mass"], psmFileReader["P-value"], psmFileReader["NIST score"], psmFileReader["Precursor Isolation m/z"], psmFileReader["Precursor Isolation Mass (Da)"], psmFileReader["Precursor Theoretical Neutral Mass (Da)"], psmFileReader["Precursor Experimental Neutral Mass (Da)"], psmFileReader["Precursor Mass Error (ppm)"], psmFileReader["Adjusted Precursor Mass Error (ppm)"], psmFileReader["Q-Value (%)"], 0,
                                    quantifiedPeptide.totalIntensity[0, NUMISOTOPES], quantifiedPeptide.totalIntensity[1, NUMISOTOPES], 0, 0, 0, 0, 
                                    0, 0, 0, 0, 0, 0, 
                                    0, 0, 0, 0, 0, 0, 0, 
                                    0, 0, 0, 0, 0, 0, "KX");
                                break;
                            case 3:
                                tagQuantWriter.WriteLine("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23},{24},{25},{26},{27},{28},{29},{30},{31},{32},{33},{34},{35},{36},{37},{38},{39},{40},{41},{42},{43},{44},{45},{46},{47},{48}",
                                    psmFileReader["Spectrum Number"], psmFileReader["Filename/id"], psmFileReader["Peptide"], psmFileReader["E-value"], psmFileReader["Mass"], psmFileReader["gi"], psmFileReader["Accession"], psmFileReader["Start"], psmFileReader["Stop"], defLineWithCommas, modsWithCommas,
                                    psmFileReader["Charge"], psmFileReader["Theo Mass"], psmFileReader["P-value"], psmFileReader["NIST score"], psmFileReader["Precursor Isolation m/z"], psmFileReader["Precursor Isolation Mass (Da)"], psmFileReader["Precursor Theoretical Neutral Mass (Da)"], psmFileReader["Precursor Experimental Neutral Mass (Da)"], psmFileReader["Precursor Mass Error (ppm)"], psmFileReader["Adjusted Precursor Mass Error (ppm)"], psmFileReader["Q-Value (%)"], 0,
                                    quantifiedPeptide.totalIntensity[0, NUMISOTOPES], quantifiedPeptide.totalIntensity[1, NUMISOTOPES], quantifiedPeptide.totalIntensity[2, NUMISOTOPES], 0, 0, 0, 
                                    0, 0, 0, 0, 0, 0, 0, 
                                    0, 0, 0, 0, 0, 0, 
                                    0, 0, 0, 0, 0, 0, "KX");
                                break;
                            case 4:
                                tagQuantWriter.WriteLine("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23},{24},{25},{26},{27},{28},{29},{30},{31},{32},{33},{34},{35},{36},{37},{38},{39},{40},{41},{42},{43},{44},{45},{46},{47},{48}",
                                    psmFileReader["Spectrum Number"], psmFileReader["Filename/id"], psmFileReader["Peptide"], psmFileReader["E-value"], psmFileReader["Mass"], psmFileReader["gi"], psmFileReader["Accession"], psmFileReader["Start"], psmFileReader["Stop"], defLineWithCommas, modsWithCommas,
                                    psmFileReader["Charge"], psmFileReader["Theo Mass"], psmFileReader["P-value"], psmFileReader["NIST score"], psmFileReader["Precursor Isolation m/z"], psmFileReader["Precursor Isolation Mass (Da)"], psmFileReader["Precursor Theoretical Neutral Mass (Da)"], psmFileReader["Precursor Experimental Neutral Mass (Da)"], psmFileReader["Precursor Mass Error (ppm)"], psmFileReader["Adjusted Precursor Mass Error (ppm)"], psmFileReader["Q-Value (%)"], 0,
                                    quantifiedPeptide.totalIntensity[0, NUMISOTOPES], quantifiedPeptide.totalIntensity[1, NUMISOTOPES], quantifiedPeptide.totalIntensity[2, NUMISOTOPES], quantifiedPeptide.totalIntensity[3, NUMISOTOPES], 0, 0, 
                                    0, 0, 0, 0, 0, 0, 0, 
                                    0, 0, 0, 0, 0, 0, 
                                    0, 0, 0, 0, 0, 0, "KX");
                                break;
                            case 6:
                                tagQuantWriter.WriteLine("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23},{24},{25},{26},{27},{28},{29},{30},{31},{32},{33},{34},{35},{36},{37},{38},{39},{40},{41},{42},{43},{44},{45},{46},{47},{48}",
                                    psmFileReader["Spectrum Number"], psmFileReader["Filename/id"], psmFileReader["Peptide"], psmFileReader["E-value"], psmFileReader["Mass"], psmFileReader["gi"], psmFileReader["Accession"], psmFileReader["Start"], psmFileReader["Stop"], defLineWithCommas, modsWithCommas,
                                    psmFileReader["Charge"], psmFileReader["Theo Mass"], psmFileReader["P-value"], psmFileReader["NIST score"], psmFileReader["Precursor Isolation m/z"], psmFileReader["Precursor Isolation Mass (Da)"], psmFileReader["Precursor Theoretical Neutral Mass (Da)"], psmFileReader["Precursor Experimental Neutral Mass (Da)"], psmFileReader["Precursor Mass Error (ppm)"], psmFileReader["Adjusted Precursor Mass Error (ppm)"], psmFileReader["Q-Value (%)"], 0,
                                    quantifiedPeptide.totalIntensity[0, NUMISOTOPES], quantifiedPeptide.totalIntensity[1, NUMISOTOPES], quantifiedPeptide.totalIntensity[2, NUMISOTOPES], quantifiedPeptide.totalIntensity[3, NUMISOTOPES], quantifiedPeptide.totalIntensity[4, NUMISOTOPES], quantifiedPeptide.totalIntensity[5, NUMISOTOPES], 
                                    0, 0, 0, 0, 0, 0, 0, 
                                    0, 0, 0, 0, 0, 0, 
                                    0, 0, 0, 0, 0, 0, "KX");
                                break;
                            default:
                                break;
                        }
                    }
                    else if (NUMCLUSTERS == 2)
                    {
                        switch (NUMISOTOPOLOGUES)
                        {
                            case 1:
                                tagQuantWriter.WriteLine("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23},{24},{25},{26},{27},{28},{29},{30},{31},{32},{33},{34},{35},{36},{37},{38},{39},{40},{41},{42},{43},{44},{45},{46},{47},{48}",
                                     psmFileReader["Spectrum Number"], psmFileReader["Filename/id"], psmFileReader["Peptide"], psmFileReader["E-value"], psmFileReader["Mass"], psmFileReader["gi"], psmFileReader["Accession"], psmFileReader["Start"], psmFileReader["Stop"], defLineWithCommas, modsWithCommas,
                                     psmFileReader["Charge"], psmFileReader["Theo Mass"], psmFileReader["P-value"], psmFileReader["NIST score"], psmFileReader["Precursor Isolation m/z"], psmFileReader["Precursor Isolation Mass (Da)"], psmFileReader["Precursor Theoretical Neutral Mass (Da)"], psmFileReader["Precursor Experimental Neutral Mass (Da)"], psmFileReader["Precursor Mass Error (ppm)"], psmFileReader["Adjusted Precursor Mass Error (ppm)"], psmFileReader["Q-Value (%)"], 0,
                                    quantifiedPeptide.totalIntensity[0, NUMISOTOPES], 0, 0, 0, 0, 0, 
                                    quantifiedPeptide.totalIntensity[1, NUMISOTOPES], 0, 0, 0, 0, 0, 0, 
                                    0, 0, 0, 0, 0, 0, 
                                    0, 0, 0, 0, 0, 0, "KX");
                                break;
                            case 2:
                                tagQuantWriter.WriteLine("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23},{24},{25},{26},{27},{28},{29},{30},{31},{32},{33},{34},{35},{36},{37},{38},{39},{40},{41},{42},{43},{44},{45},{46},{47},{48}",
                                    psmFileReader["Spectrum Number"], psmFileReader["Filename/id"], psmFileReader["Peptide"], psmFileReader["E-value"], psmFileReader["Mass"], psmFileReader["gi"], psmFileReader["Accession"], psmFileReader["Start"], psmFileReader["Stop"], defLineWithCommas, modsWithCommas,
                                    psmFileReader["Charge"], psmFileReader["Theo Mass"], psmFileReader["P-value"], psmFileReader["NIST score"], psmFileReader["Precursor Isolation m/z"], psmFileReader["Precursor Isolation Mass (Da)"], psmFileReader["Precursor Theoretical Neutral Mass (Da)"], psmFileReader["Precursor Experimental Neutral Mass (Da)"], psmFileReader["Precursor Mass Error (ppm)"], psmFileReader["Adjusted Precursor Mass Error (ppm)"], psmFileReader["Q-Value (%)"], 0,
                                    quantifiedPeptide.totalIntensity[0, NUMISOTOPES], quantifiedPeptide.totalIntensity[1, NUMISOTOPES], 0, 0, 0, 0, 
                                    quantifiedPeptide.totalIntensity[2, NUMISOTOPES], quantifiedPeptide.totalIntensity[3, NUMISOTOPES], 0, 0, 0, 0, 0, 
                                    0, 0, 0, 0, 0, 0, 
                                    0, 0, 0, 0, 0, 0, "KX");
                                break;
                            case 3:
                                tagQuantWriter.WriteLine("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23},{24},{25},{26},{27},{28},{29},{30},{31},{32},{33},{34},{35},{36},{37},{38},{39},{40},{41},{42},{43},{44},{45},{46},{47},{48}",
                                    psmFileReader["Spectrum Number"], psmFileReader["Filename/id"], psmFileReader["Peptide"], psmFileReader["E-value"], psmFileReader["Mass"], psmFileReader["gi"], psmFileReader["Accession"], psmFileReader["Start"], psmFileReader["Stop"], defLineWithCommas, modsWithCommas,
                                    psmFileReader["Charge"], psmFileReader["Theo Mass"], psmFileReader["P-value"], psmFileReader["NIST score"], psmFileReader["Precursor Isolation m/z"], psmFileReader["Precursor Isolation Mass (Da)"], psmFileReader["Precursor Theoretical Neutral Mass (Da)"], psmFileReader["Precursor Experimental Neutral Mass (Da)"], psmFileReader["Precursor Mass Error (ppm)"], psmFileReader["Adjusted Precursor Mass Error (ppm)"], psmFileReader["Q-Value (%)"], 0,
                                    quantifiedPeptide.totalIntensity[0, NUMISOTOPES], quantifiedPeptide.totalIntensity[1, NUMISOTOPES], quantifiedPeptide.totalIntensity[2, NUMISOTOPES], 0, 0, 0, 
                                    quantifiedPeptide.totalIntensity[3, NUMISOTOPES], quantifiedPeptide.totalIntensity[4, NUMISOTOPES], quantifiedPeptide.totalIntensity[5, NUMISOTOPES], 0, 0, 0, 0, 
                                    0, 0, 0, 0, 0, 0, 
                                    0, 0, 0, 0, 0, 0, "KX");
                                break;
                            case 4:
                                tagQuantWriter.WriteLine("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23},{24},{25},{26},{27},{28},{29},{30},{31},{32},{33},{34},{35},{36},{37},{38},{39},{40},{41},{42},{43},{44},{45},{46},{47},{48}",
                                    psmFileReader["Spectrum Number"], psmFileReader["Filename/id"], psmFileReader["Peptide"], psmFileReader["E-value"], psmFileReader["Mass"], psmFileReader["gi"], psmFileReader["Accession"], psmFileReader["Start"], psmFileReader["Stop"], defLineWithCommas, modsWithCommas,
                                    psmFileReader["Charge"], psmFileReader["Theo Mass"], psmFileReader["P-value"], psmFileReader["NIST score"], psmFileReader["Precursor Isolation m/z"], psmFileReader["Precursor Isolation Mass (Da)"], psmFileReader["Precursor Theoretical Neutral Mass (Da)"], psmFileReader["Precursor Experimental Neutral Mass (Da)"], psmFileReader["Precursor Mass Error (ppm)"], psmFileReader["Adjusted Precursor Mass Error (ppm)"], psmFileReader["Q-Value (%)"], 0,
                                    quantifiedPeptide.totalIntensity[0, NUMISOTOPES], quantifiedPeptide.totalIntensity[1, NUMISOTOPES], quantifiedPeptide.totalIntensity[2, NUMISOTOPES], quantifiedPeptide.totalIntensity[3, NUMISOTOPES], 0, 0, 
                                    quantifiedPeptide.totalIntensity[4, NUMISOTOPES], quantifiedPeptide.totalIntensity[5, NUMISOTOPES], quantifiedPeptide.totalIntensity[6, NUMISOTOPES], quantifiedPeptide.totalIntensity[7, NUMISOTOPES], 0, 0, 0, 
                                    0, 0, 0, 0, 0, 0, 
                                    0, 0, 0, 0, 0, 0, "KX");
                                break;
                            case 6:
                                tagQuantWriter.WriteLine("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23},{24},{25},{26},{27},{28},{29},{30},{31},{32},{33},{34},{35},{36},{37},{38},{39},{40},{41},{42},{43},{44},{45},{46},{47},{48}",
                                    psmFileReader["Spectrum Number"], psmFileReader["Filename/id"], psmFileReader["Peptide"], psmFileReader["E-value"], psmFileReader["Mass"], psmFileReader["gi"], psmFileReader["Accession"], psmFileReader["Start"], psmFileReader["Stop"], defLineWithCommas, modsWithCommas,
                                    psmFileReader["Charge"], psmFileReader["Theo Mass"], psmFileReader["P-value"], psmFileReader["NIST score"], psmFileReader["Precursor Isolation m/z"], psmFileReader["Precursor Isolation Mass (Da)"], psmFileReader["Precursor Theoretical Neutral Mass (Da)"], psmFileReader["Precursor Experimental Neutral Mass (Da)"], psmFileReader["Precursor Mass Error (ppm)"], psmFileReader["Adjusted Precursor Mass Error (ppm)"], psmFileReader["Q-Value (%)"], 0,
                                    quantifiedPeptide.totalIntensity[0, NUMISOTOPES], quantifiedPeptide.totalIntensity[1, NUMISOTOPES], quantifiedPeptide.totalIntensity[2, NUMISOTOPES], quantifiedPeptide.totalIntensity[3, NUMISOTOPES], quantifiedPeptide.totalIntensity[4, NUMISOTOPES], quantifiedPeptide.totalIntensity[5, NUMISOTOPES], 
                                    quantifiedPeptide.totalIntensity[6, NUMISOTOPES], quantifiedPeptide.totalIntensity[7, NUMISOTOPES], quantifiedPeptide.totalIntensity[8, NUMISOTOPES], quantifiedPeptide.totalIntensity[9, NUMISOTOPES], quantifiedPeptide.totalIntensity[10, NUMISOTOPES], quantifiedPeptide.totalIntensity[11, NUMISOTOPES], 0, 
                                    0, 0, 0, 0, 0, 0, 
                                    0, 0, 0, 0, 0, 0, "KX");
                                break;
                            default:
                                break;
                        }
                    }
                    else if (NUMCLUSTERS == 3)
                    {
                        switch (NUMISOTOPOLOGUES)
                        {
                            case 1:
                                tagQuantWriter.WriteLine("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23},{24},{25},{26},{27},{28},{29},{30},{31},{32},{33},{34},{35},{36},{37},{38},{39},{40},{41},{42},{43},{44},{45},{46},{47},{48}",
                                    psmFileReader["Spectrum Number"], psmFileReader["Filename/id"], psmFileReader["Peptide"], psmFileReader["E-value"], psmFileReader["Mass"], psmFileReader["gi"], psmFileReader["Accession"], psmFileReader["Start"], psmFileReader["Stop"], defLineWithCommas, modsWithCommas,
                                    psmFileReader["Charge"], psmFileReader["Theo Mass"], psmFileReader["P-value"], psmFileReader["NIST score"], psmFileReader["Precursor Isolation m/z"], psmFileReader["Precursor Isolation Mass (Da)"], psmFileReader["Precursor Theoretical Neutral Mass (Da)"], psmFileReader["Precursor Experimental Neutral Mass (Da)"], psmFileReader["Precursor Mass Error (ppm)"], psmFileReader["Adjusted Precursor Mass Error (ppm)"], psmFileReader["Q-Value (%)"], 0,
                                    quantifiedPeptide.totalIntensity[0, NUMISOTOPES], 0, 0, 0, 0, 0,
                                    quantifiedPeptide.totalIntensity[1, NUMISOTOPES], 0, 0, 0, 0, 0, 0,
                                    quantifiedPeptide.totalIntensity[2, NUMISOTOPES], 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, "KX");
                                break;
                            case 2:
                                tagQuantWriter.WriteLine("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23},{24},{25},{26},{27},{28},{29},{30},{31},{32},{33},{34},{35},{36},{37},{38},{39},{40},{41},{42},{43},{44},{45},{46},{47},{48}",
                                    psmFileReader["Spectrum Number"], psmFileReader["Filename/id"], psmFileReader["Peptide"], psmFileReader["E-value"], psmFileReader["Mass"], psmFileReader["gi"], psmFileReader["Accession"], psmFileReader["Start"], psmFileReader["Stop"], defLineWithCommas, modsWithCommas,
                                    psmFileReader["Charge"], psmFileReader["Theo Mass"], psmFileReader["P-value"], psmFileReader["NIST score"], psmFileReader["Precursor Isolation m/z"], psmFileReader["Precursor Isolation Mass (Da)"], psmFileReader["Precursor Theoretical Neutral Mass (Da)"], psmFileReader["Precursor Experimental Neutral Mass (Da)"], psmFileReader["Precursor Mass Error (ppm)"], psmFileReader["Adjusted Precursor Mass Error (ppm)"], psmFileReader["Q-Value (%)"], 0,
                                    quantifiedPeptide.totalIntensity[0, NUMISOTOPES], quantifiedPeptide.totalIntensity[1, NUMISOTOPES], 0, 0, 0, 0,
                                    quantifiedPeptide.totalIntensity[2, NUMISOTOPES], quantifiedPeptide.totalIntensity[3, NUMISOTOPES], 0, 0, 0, 0, 0,
                                    quantifiedPeptide.totalIntensity[4, NUMISOTOPES], quantifiedPeptide.totalIntensity[5, NUMISOTOPES], 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, "KX");
                                break;
                            case 3:
                                tagQuantWriter.WriteLine("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23},{24},{25},{26},{27},{28},{29},{30},{31},{32},{33},{34},{35},{36},{37},{38},{39},{40},{41},{42},{43},{44},{45},{46},{47},{48}",
                                    psmFileReader["Spectrum Number"], psmFileReader["Filename/id"], psmFileReader["Peptide"], psmFileReader["E-value"], psmFileReader["Mass"], psmFileReader["gi"], psmFileReader["Accession"], psmFileReader["Start"], psmFileReader["Stop"], defLineWithCommas, modsWithCommas,
                                    psmFileReader["Charge"], psmFileReader["Theo Mass"], psmFileReader["P-value"], psmFileReader["NIST score"], psmFileReader["Precursor Isolation m/z"], psmFileReader["Precursor Isolation Mass (Da)"], psmFileReader["Precursor Theoretical Neutral Mass (Da)"], psmFileReader["Precursor Experimental Neutral Mass (Da)"], psmFileReader["Precursor Mass Error (ppm)"], psmFileReader["Adjusted Precursor Mass Error (ppm)"], psmFileReader["Q-Value (%)"], 0,
                                    quantifiedPeptide.totalIntensity[0, NUMISOTOPES], quantifiedPeptide.totalIntensity[1, NUMISOTOPES], quantifiedPeptide.totalIntensity[2, NUMISOTOPES], 0, 0, 0,
                                    quantifiedPeptide.totalIntensity[3, NUMISOTOPES], quantifiedPeptide.totalIntensity[4, NUMISOTOPES], quantifiedPeptide.totalIntensity[5, NUMISOTOPES], 0, 0, 0, 0,
                                    quantifiedPeptide.totalIntensity[6, NUMISOTOPES], quantifiedPeptide.totalIntensity[7, NUMISOTOPES], quantifiedPeptide.totalIntensity[8, NUMISOTOPES], 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, "KX");
                                break;
                            case 4:
                                tagQuantWriter.WriteLine("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23},{24},{25},{26},{27},{28},{29},{30},{31},{32},{33},{34},{35},{36},{37},{38},{39},{40},{41},{42},{43},{44},{45},{46},{47},{48}",
                                    psmFileReader["Spectrum Number"], psmFileReader["Filename/id"], psmFileReader["Peptide"], psmFileReader["E-value"], psmFileReader["Mass"], psmFileReader["gi"], psmFileReader["Accession"], psmFileReader["Start"], psmFileReader["Stop"], defLineWithCommas, modsWithCommas,
                                    psmFileReader["Charge"], psmFileReader["Theo Mass"], psmFileReader["P-value"], psmFileReader["NIST score"], psmFileReader["Precursor Isolation m/z"], psmFileReader["Precursor Isolation Mass (Da)"], psmFileReader["Precursor Theoretical Neutral Mass (Da)"], psmFileReader["Precursor Experimental Neutral Mass (Da)"], psmFileReader["Precursor Mass Error (ppm)"], psmFileReader["Adjusted Precursor Mass Error (ppm)"], psmFileReader["Q-Value (%)"], 0,
                                    quantifiedPeptide.totalIntensity[0, NUMISOTOPES], quantifiedPeptide.totalIntensity[1, NUMISOTOPES], quantifiedPeptide.totalIntensity[2, NUMISOTOPES], quantifiedPeptide.totalIntensity[3, NUMISOTOPES], 0, 0,
                                    quantifiedPeptide.totalIntensity[4, NUMISOTOPES], quantifiedPeptide.totalIntensity[5, NUMISOTOPES], quantifiedPeptide.totalIntensity[6, NUMISOTOPES], quantifiedPeptide.totalIntensity[7, NUMISOTOPES], 0, 0, 0,
                                    quantifiedPeptide.totalIntensity[8, NUMISOTOPES], quantifiedPeptide.totalIntensity[9, NUMISOTOPES], quantifiedPeptide.totalIntensity[10, NUMISOTOPES], quantifiedPeptide.totalIntensity[11, NUMISOTOPES], 0, 0,
                                    0, 0, 0, 0, 0, 0, "KX");
                                break;
                            case 6:
                                tagQuantWriter.WriteLine("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23},{24},{25},{26},{27},{28},{29},{30},{31},{32},{33},{34},{35},{36},{37},{38},{39},{40},{41},{42},{43},{44},{45},{46},{47},{48}",
                                    psmFileReader["Spectrum Number"], psmFileReader["Filename/id"], psmFileReader["Peptide"], psmFileReader["E-value"], psmFileReader["Mass"], psmFileReader["gi"], psmFileReader["Accession"], psmFileReader["Start"], psmFileReader["Stop"], defLineWithCommas, modsWithCommas,
                                    psmFileReader["Charge"], psmFileReader["Theo Mass"], psmFileReader["P-value"], psmFileReader["NIST score"], psmFileReader["Precursor Isolation m/z"], psmFileReader["Precursor Isolation Mass (Da)"], psmFileReader["Precursor Theoretical Neutral Mass (Da)"], psmFileReader["Precursor Experimental Neutral Mass (Da)"], psmFileReader["Precursor Mass Error (ppm)"], psmFileReader["Adjusted Precursor Mass Error (ppm)"], psmFileReader["Q-Value (%)"], 0,
                                    quantifiedPeptide.totalIntensity[0, NUMISOTOPES], quantifiedPeptide.totalIntensity[1, NUMISOTOPES], quantifiedPeptide.totalIntensity[2, NUMISOTOPES], quantifiedPeptide.totalIntensity[3, NUMISOTOPES], quantifiedPeptide.totalIntensity[4, NUMISOTOPES], quantifiedPeptide.totalIntensity[5, NUMISOTOPES],
                                    quantifiedPeptide.totalIntensity[6, NUMISOTOPES], quantifiedPeptide.totalIntensity[7, NUMISOTOPES], quantifiedPeptide.totalIntensity[8, NUMISOTOPES], quantifiedPeptide.totalIntensity[9, NUMISOTOPES], quantifiedPeptide.totalIntensity[10, NUMISOTOPES], quantifiedPeptide.totalIntensity[11, NUMISOTOPES], 0,
                                    quantifiedPeptide.totalIntensity[12, NUMISOTOPES], quantifiedPeptide.totalIntensity[13, NUMISOTOPES], quantifiedPeptide.totalIntensity[14, NUMISOTOPES], quantifiedPeptide.totalIntensity[15, NUMISOTOPES], quantifiedPeptide.totalIntensity[16, NUMISOTOPES], quantifiedPeptide.totalIntensity[17, NUMISOTOPES],
                                    0, 0, 0, 0, 0, 0, "KX");
                                break;
                            default:
                                break;
                        }
                    }
                }
            }
            tagQuantWriter.Close();
        }

        /* Reads in peptide information from a CSV file
         * Necessary information: spectrum number, charge, peptide sequence, E-value, raw file
         */
        private void readCsvInputFile(Dictionary<string, PeptideID> allPeptides, string file)
        {
            HashSet<string> rawFiles = new HashSet<string>();

            //Cycle through .csv file to make a list of identified peptides and properties
            CsvReader psmFileReader = new CsvReader(new StreamReader(file), true);
            CsvReader injectionTimesReader;

            using (psmFileReader)
            {
                PeptideID newPeptide;
                MSDataFile rawFile = null;
                WriteMessage("uploading peptide IDs");
                while (psmFileReader.ReadNextRecord())
                {
                    string basePathName = psmFileReader["Filename/id"].Substring(0, psmFileReader["Filename/id"].IndexOf("."));

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

                    newPeptide = new PeptideID(int.Parse(psmFileReader["Spectrum number"]),
                                int.Parse(psmFileReader["Charge"]), double.Parse(psmFileReader["E-value"]),
                                psmFileReader["Peptide"], rawFile, psmFileReader["Mods"]);
                    checkAdd(allPeptides, newPeptide, (MSDataFile)rawFile, double.Parse(psmFileReader["E-value"]));                 
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
        private double calculateCoalescenceThreshold(string file, Dictionary<string, PeptideID> allPeptides, double binSize, List<CoalescenceCheck> list)
        {
            double intensityCutOff = 10000000.0;
            double intensityBinSize = binSize;
            CoalescenceCheck check;

            StreamWriter coalWriter = new StreamWriter(Path.Combine(outputFolderBox.Text, Path.GetFileNameWithoutExtension(file) + "_coalescencePlot.csv"));
            string header = ("Intensity Bin, # Peptides, Average Missing Channel Frequency");
            coalWriter.WriteLine(header);
            WriteMessage("writing output");
    
            //Retrieve each peptide's maximum intensity & missing channel frequency
            foreach (PeptideID uniquePeptide in allPeptides.Values)
            {
                double maximumIncompleteIntensity;
                double maximumCompleteIntensity;
                double maximumIntensity;
                double missingChannelFrequency;
                if (NUMISOTOPOLOGUES > 1)
                {
                    for (int c = 0; c < NUMCLUSTERS; c++)
                    {
                        maximumIncompleteIntensity = uniquePeptide.log10MaxIntensity[c];
                        maximumCompleteIntensity = uniquePeptide.log10MaxCompleteIntensity[c];
                        if (maximumCompleteIntensity >= maximumIncompleteIntensity) maximumIntensity = maximumCompleteIntensity;
                        else
                        {
                            maximumIntensity = maximumIncompleteIntensity;
                        }
                        missingChannelFrequency = uniquePeptide.missingChannelFrequency[c];
                        if (maximumIntensity > 0 && missingChannelFrequency >= 0)
                        {
                            check = new CoalescenceCheck(maximumIntensity, missingChannelFrequency);
                            list.Add(check);
                        }
                    }
                }
                else
                {

                }
            }

            List<CoalescenceCheck> sortedList = list.OrderBy(coalCheck => coalCheck.intensity).ToList();

            // Set up the intensity bins (log10)
            int intensityCount = sortedList.Count;
            double minIntensity = sortedList[0].intensity;
            double maxIntensity = sortedList[intensityCount - 1].intensity;
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

                coalWriter.WriteLine("{0},{1},{2}", intensityValue, count, averageFrequency);
                
                // Check if bin average frequency is one of the 5 lowest out of all intensity bins
                if (averageFrequency < sortedLowestFreqs.ElementAt(4).Key && averageFrequency > 0)
                {
                    intensityValue = (double)(i * intensityBinSize) + firstBinMaximum;
                    sortedLowestFreqs.Add(averageFrequency, intensityValue);
                    sortedLowestFreqs.RemoveAt(5);
                }
                binAverageFrequency[i] = averageFrequency;
            }
            coalWriter.Close();

            // Calculate and return the average intensity producing the minimum missing channel frequency
            sum = 0;
            foreach (KeyValuePair<double, double> pair in sortedLowestFreqs)
            {
                sum += pair.Value;
            }

            intensityCutOff = System.Math.Pow(10.0, sum / 5.0);

            coalWriter.Close();
            
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
            OnMessage += Form1_OnMessage;
            Thread thread = new Thread(Run);
            thread.IsBackground = true;
            thread.Start();
            start.Enabled = true;            
        }

        private delegate void OnMessageDelegate(string msg);

        private void Form1_OnMessage(object sender, MessageEventArgs e)
        {
            WriteMessage(e.Message);            
        }

        private void Run()
        {
            foreach (string file in listBox1.Items)
            {
                Run(file);
            }
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

                //int index = files[0].LastIndexOf("\\");
                //rawFileBox.Text = files[0].Substring(0, index + 1);
                //outputFolderBox.Text = rawFileBox.Text;
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
