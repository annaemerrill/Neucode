using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Coon;
using CoonThermo.IO;
using Coon.Spectra;


namespace HILAC
{
    /* OMSSA-identified or targeted peptide to be quantified by HILAC or SILAC
     * Keeps track of a list of HILAC pairs for quantitation
    */
    class AveragePeptideID
    {
        //mTRAQ modifications & constants
        public static Mass lightmTRAQ = new Mass(140.0953, 140.0953, 140.0953);
        public static Mass lightmTRAQK = new Mass(148.1275, 148.1275, 148.1275);
        public static Modification lightmTRAQNTerm = new Modification("mTRAQ Light N-term", lightmTRAQ);
        public static Modification lightmTRAQLysine = new Modification("mTRAQ Light Lysine", lightmTRAQK);
        //public static Modification lightmTRAQTyrosine = new Modification("mTRAQ Light Tyrosine", lightmTRAQ);
        public static Mass mediummTRAQ = new Mass(144.1024, 144.1024, 144.1024);
        public static Mass mediummTRAQK = new Mass(152.1346, 152.1346, 152.1346);
        public static Modification mediummTRAQNTerm = new Modification("mTRAQ Medium N-term", mediummTRAQ);
        public static Modification mediummTRAQLysine = new Modification("mTRAQ Medium Lysine", mediummTRAQK);
        //public static Modification mediummTRAQTyrosine = new Modification("mTRAQ Medium Tyrosine", mediummTRAQ);
        public static Mass heavymTRAQ = new Mass(148.1040, 148.1040, 148.1040);
        public static Mass heavymTRAQK = new Mass(156.1362, 156.1362, 156.1362);
        public static Modification heavymTRAQNTerm = new Modification("mTRAQ Heavy N-term", heavymTRAQ);
        public static Modification heavymTRAQLysine = new Modification("mTRAQ Heavy Lysine", heavymTRAQK);
        //public static Modification heavymTRAQTyrosine = new Modification("mTRAQ Heavy Tyrosine", heavymTRAQ);
        public static ChemicalModification CAM = new ChemicalModification("C2H3NO");
        public static ChemicalModification METOX = new ChemicalModification("O");

        //Used only for 6-plex HILAC parent
        public AveragePeptideID parent; 
        public bool light; 
        public bool medium;
        public bool heavy;
        public Peptide peptide;
        public AveragePeptideID mTRAQLight; //mTRAQ Light
        public AveragePeptideID mTRAQMedium; //mTRAQ Medium
        public AveragePeptideID mTRAQHeavy; //mTRAQ Heavy
        public bool identified;
        
        //Used for all HILAC and SILAC
        public SortedList<double, PeptideSpectralMatch> PSMs;
        public string sequence;
        public int charge
        {
            get
            {
                if (PSMs != null && PSMs.Count > 0)
                {
                    return PSMs.ElementAt(0).Value.Charge;
                }
                return 0;
            }
        }
        public List<int> chargeStates;
        public double theoreticalMass;
        public double experimentalMass;
        public double lightTotalIntensity;
        public double heavyTotalIntensity;
        public List<HILACPair> hILACPairs;
        public List<HILACPair> alignedhILACPairs;
        public List<HILACPair> filteredhILACPairs;
        public SortedList<double, HILACPairIsotope> quantifiableHILACPairIsotopes;
        public MetabolicForm maxLight;
        public MetabolicForm maxHeavy;
        public double lightIntensityCutOff
        {
            get
            {
                if (maxLight != null)
                {
                    return (Form1.QUANTINTENSITY * maxLight.intensity);
                }
                return double.NaN;
            }
        }
        public double heavyIntensityCutOff
        {
            get
            {
                if (maxHeavy != null)
                {
                    return (Form1.QUANTINTENSITY * maxHeavy.intensity);
                }
                return double.NaN;
            }
        }
        public int peakMaxShift //# of MS1 scans between heavy peak max and light peak max ((-) --> heavy elutes before light; (+) --> heavy elutes after light)
        {
            get
            {
                if (maxLight != null && maxHeavy != null)
                {
                    if (maxHeavy.pair.noiseBandCapped || maxLight.pair.noiseBandCapped)
                    {
                        return 0;
                    }
                    else
                    {
                        return (hILACPairs.IndexOf(maxHeavy.pair.pair) - hILACPairs.IndexOf(maxLight.pair.pair));
                    }
                }
                return 0;
            }
            set { }
        }
        public int firstScanNumber
        {
            get
            {
                if (PSMs != null && PSMs.Count > 0)
                {
                    int number = PSMs.ElementAt(0).Value.ScanNumber;
                    for (int i = 1; i < PSMs.Count(); i++)
                    {
                        if (PSMs.ElementAt(i).Value.ScanNumber < number)
                        {
                            number = PSMs.ElementAt(i).Value.ScanNumber;
                        }
                    }
                    return number;
                }
                return 0;
            }
        }
        public int lastScanNumber
        {
            get
            {
                if (PSMs != null && PSMs.Count > 0)
                {
                    int number = PSMs.ElementAt(0).Value.ScanNumber;
                    for (int i = 1; i < PSMs.Count(); i++)
                    {
                        if (PSMs.ElementAt(i).Value.ScanNumber > number)
                        {
                            number = PSMs.ElementAt(i).Value.ScanNumber;
                        }
                    }
                    return number;
                }
                return 0;
            }
        }

        public int lysineCount;
        public double calculatedIntraMZDifference;

        public Range rTRange;
        
        public Range precursorRange0;
        public Range precursorRange1;
        public Range precursorRange2;
        public Range precursorRange3;
        public Range fullScanRange;
        public double medianHeavyToLightRatio;
        public double slopeHeavyToLightRatio;
        /*public ThermoRawFileScan MS
        {
            get
            {
                if (MSMS != null)
                {
                    return MSMS.GetPrecursorRawScan();
                }
                return null;
            }
        }*/
        public bool noiseBandCapped;
        public bool reQuantify;

        //public ThermoRawFileScan MSMS;
        public int reasonNQ;

        //Constructor for duplex HILAC or SILAC
        public AveragePeptideID(int scanNumber, int charge, double theoMass, double expMass, double eValue, string sequence, double lysineVersion)
        {
            PSMs = new SortedList<double, PeptideSpectralMatch>();
            PSMs.Add(eValue, new PeptideSpectralMatch(scanNumber, charge, sequence, Form1.RAWFILE[scanNumber].Spectrum));
            theoreticalMass = theoMass;
            experimentalMass = expMass;
            this.sequence = sequence;
            lightTotalIntensity = 0;
            heavyTotalIntensity = 0;
            hILACPairs = new List<HILACPair>();
            countLysines();
            maxLight = null;
            maxHeavy = null;
            noiseBandCapped = false;
            reQuantify = false;
        }

        //Constructor for 6-plex HILAC
        public AveragePeptideID(int scanNumber, int charge, double theoMass, double expMass, double eValue, string sequence, string modifications)
        {
            PSMs = new SortedList<double, PeptideSpectralMatch>();
            light = false;
            medium = false;
            heavy = false;
            theoreticalMass = theoMass;
            experimentalMass = expMass;
            this.sequence = sequence;
            countLysines();
            lightTotalIntensity = 0;
            heavyTotalIntensity = 0;           

            Peptide peptide = new Peptide(sequence, modifications); //Figure out which version was identified
            peptide.SetFixedModifications('C', CAM);

            Peptide checkLight = new Peptide(peptide);
            Peptide checkMedium = new Peptide(peptide);
            Peptide checkHeavy = new Peptide(peptide);

            //Try setting identified peptide as "light", "medium", or "heavy" & check error from theoretical mass
            checkLight.SetFixedModifications('K', lightmTRAQLysine);
            checkLight.SetFixedModification(1, lightmTRAQNTerm);
            
            checkMedium.SetFixedModifications('K', mediummTRAQLysine);
            checkMedium.SetFixedModification(1, mediummTRAQNTerm);

            checkHeavy.SetFixedModifications('K', heavymTRAQLysine);
            checkHeavy.SetFixedModification(1, heavymTRAQNTerm);

            double checkLightError = Tolerance.GetError(theoreticalMass, checkLight.Mass.MonoisotopicMass, ToleranceType.PPM);
            double checkMediumError = Tolerance.GetError(theoreticalMass, checkMedium.Mass.MonoisotopicMass, ToleranceType.PPM);
            double checkHeavyError = Tolerance.GetError(theoreticalMass, checkHeavy.Mass.MonoisotopicMass, ToleranceType.PPM);

            //Set as mTRAQ version with smallest mass error
            if (Math.Abs(checkLightError) <= Math.Abs(checkMediumError) && Math.Abs(checkLightError) <= Math.Abs(checkHeavyError)) //Identification is "light"
            {
                light = true;
                mTRAQLight = new AveragePeptideID(scanNumber, charge, checkLight.Mass.MonoisotopicMass, experimentalMass, eValue, sequence, Form1.HILACDIFF);
                mTRAQLight.parent = this;
                mTRAQLight.identified = true;
                PSMs.Add(eValue, new PeptideSpectralMatch(scanNumber, charge, checkLight, Form1.RAWFILE[scanNumber].Spectrum));
                mTRAQMedium = new AveragePeptideID(scanNumber, charge, checkMedium.Mass.MonoisotopicMass, experimentalMass, eValue, sequence, Form1.HILACDIFF);
                mTRAQMedium.parent = this;
                mTRAQHeavy = new AveragePeptideID(scanNumber, charge, checkHeavy.Mass.MonoisotopicMass, experimentalMass, eValue, sequence, Form1.HILACDIFF);
                mTRAQHeavy.parent = this;
            }
            else if (Math.Abs(checkMediumError) <= Math.Abs(checkHeavyError) && Math.Abs(checkMediumError) <= Math.Abs(checkLightError)) //Identification is "medium"
            {
                medium = true;
                mTRAQMedium = new AveragePeptideID(scanNumber, charge, checkMedium.Mass.MonoisotopicMass, experimentalMass, eValue, sequence, Form1.HILACDIFF);
                mTRAQMedium.parent = this;
                mTRAQMedium.identified = true;
                PSMs.Add(eValue, new PeptideSpectralMatch(scanNumber, charge, checkMedium, Form1.RAWFILE[scanNumber].Spectrum));
                mTRAQLight = new AveragePeptideID(scanNumber, charge, checkLight.Mass.MonoisotopicMass, experimentalMass, eValue, sequence, Form1.HILACDIFF);
                mTRAQLight.parent = this;
                mTRAQHeavy = new AveragePeptideID(scanNumber, charge, checkHeavy.Mass.MonoisotopicMass, experimentalMass, eValue, sequence, Form1.HILACDIFF);
                mTRAQHeavy.parent = this;
            }
            else //Identification is "heavy"
            {
                heavy = true;
                mTRAQHeavy = new AveragePeptideID(scanNumber, charge, checkHeavy.Mass.MonoisotopicMass, experimentalMass, eValue, sequence, Form1.HILACDIFF);
                mTRAQHeavy.parent = this;
                mTRAQHeavy.identified = true;
                PSMs.Add(eValue, new PeptideSpectralMatch(scanNumber, charge, checkHeavy, Form1.RAWFILE[scanNumber].Spectrum));
                mTRAQLight = new AveragePeptideID(scanNumber, charge, checkLight.Mass.MonoisotopicMass, experimentalMass, eValue, sequence, Form1.HILACDIFF);
                mTRAQLight.parent = this;
                mTRAQMedium = new AveragePeptideID(scanNumber, charge, checkMedium.Mass.MonoisotopicMass, experimentalMass, eValue, sequence, Form1.HILACDIFF);
                mTRAQMedium.parent = this;
            }        
        }

        public void countLysines()
        {
            lysineCount = 0;
            string peptid = sequence;
            char[] residues = peptid.ToCharArray();
            for (int i = 0; i < residues.Length; i++)
            {
                if (residues[i] == 'K')
                {
                    lysineCount++;
                }
            }
        }

        //Calculate the MS1 scan range to use for quantification (include all PSMs; extending above and below according to the entered RT window)
        public void calculateScanRange()
        {
            rTRange = new Range(Form1.RAWFILE[firstScanNumber].GetPrecursorRawScan().ScanTime - Form1.RTWINDOW, Form1.RAWFILE[lastScanNumber].GetPrecursorRawScan().ScanTime + Form1.RTWINDOW);

            double firstScanNum;
            double lastScanNum;
            double minValue;
            double maxValue;
            ThermoRawFileScan current = Form1.RAWFILE[firstScanNumber].GetPrecursorRawScan();
            double firstTime = Form1.RAWFILE[Form1.RAWFILE.FirstScan].ScanTime;
            double lastTime = Form1.RAWFILE[Form1.RAWFILE.LastScan].ScanTime;

            if (rTRange.MinValue < firstTime)
            {
                minValue = firstTime;
            }
            else
            {
                minValue = rTRange.MinValue;
            }
            if (rTRange.MaxValue > lastTime)
            {
                maxValue = lastTime;
            }
            else
            {
                maxValue = rTRange.MaxValue;
            }

            if (Form1.HILAC)
            {
                while (current.GetPreviousMsnScan(1).GetPreviousMsnScan(1) != null && current.GetPreviousMsnScan(1).GetPreviousMsnScan(1).ScanTime > rTRange.MinValue)
                {
                    current = current.GetPreviousMsnScan(1).GetPreviousMsnScan(1);
                }
                firstScanNum = current.ScanNum;

                current = Form1.RAWFILE[lastScanNumber].GetPrecursorRawScan();
                while (current.GetNextMsnScan(1).GetNextMsnScan(1) != null && current.GetNextMsnScan(1).GetNextMsnScan(1).ScanTime < rTRange.MaxValue)
                {
                    current = current.GetNextMsnScan(1).GetNextMsnScan(1);
                }
                lastScanNum = current.ScanNum;
            }
            else
            {
                while (current.GetPreviousMsnScan(1) != null && current.GetPreviousMsnScan(1).ScanTime > minValue)
                {
                    current = current.GetPreviousMsnScan(1);
                }
                firstScanNum = current.ScanNum;

                current = Form1.RAWFILE[lastScanNumber].GetPrecursorRawScan();
                while (current.GetNextMsnScan(1) != null && current.GetNextMsnScan(1).ScanTime < maxValue)
                {
                    current = current.GetNextMsnScan(1);
                }
                lastScanNum = current.ScanNum;
            }
            fullScanRange = new Range(firstScanNum, lastScanNum);
        }

        //Compile a list of all charge states that led to an identification
        public void findChargeStates()
        {
            chargeStates = new List<int>();

            for (int i = 0; i < PSMs.Count; i++)
            {
                if (!chargeStates.Contains(PSMs.Values.ElementAt(i).Charge))
                {
                    chargeStates.Add(PSMs.Values.ElementAt(i).Charge);
                }
            }
        }

        public void calculatePrecursorRanges()
        {
            double lysineVersion;

            if (Form1.HILAC)
            {
                lysineVersion = 0.036;
            }
            else
            {
                lysineVersion = 8.0142;
            }

            if (!Form1.ULTRAHILAC)
            {
                calculatedIntraMZDifference = ((double)lysineCount * lysineVersion) / (double)charge;
                precursorRange0 = new Range(Mass.MZFromMass(experimentalMass, charge), new Tolerance(Form1.SEARCHPPM, ToleranceType.PPM));
                precursorRange1 = new Range(Mass.MZFromMass(experimentalMass + Constants.Neutron, charge), new Tolerance(Form1.SEARCHPPM, ToleranceType.PPM));
                precursorRange2 = new Range(Mass.MZFromMass(experimentalMass + (2.0 * Constants.Neutron), charge), new Tolerance(Form1.SEARCHPPM, ToleranceType.PPM));
                precursorRange3 = new Range(Mass.MZFromMass(experimentalMass + (3.0 * Constants.Neutron), charge), new Tolerance(Form1.SEARCHPPM, ToleranceType.PPM));
            }
            else
            {
                mTRAQLight.calculatedIntraMZDifference = ((double)lysineCount * lysineVersion) / (double)charge;
                mTRAQLight.precursorRange0 = new Range(Mass.MZFromMass(mTRAQLight.experimentalMass, charge), new Tolerance(Form1.SEARCHPPM, ToleranceType.PPM));
                mTRAQLight.precursorRange1 = new Range(Mass.MZFromMass(mTRAQLight.experimentalMass + Constants.Neutron, charge), new Tolerance(Form1.SEARCHPPM, ToleranceType.PPM));
                mTRAQLight.precursorRange2 = new Range(Mass.MZFromMass(mTRAQLight.experimentalMass + (2.0 * Constants.Neutron), charge), new Tolerance(Form1.SEARCHPPM, ToleranceType.PPM));
                mTRAQLight.precursorRange3 = new Range(Mass.MZFromMass(mTRAQLight.experimentalMass + (3.0 * Constants.Neutron), charge), new Tolerance(Form1.SEARCHPPM, ToleranceType.PPM));

                mTRAQMedium.calculatedIntraMZDifference = ((double)lysineCount * lysineVersion) / (double)charge;
                mTRAQMedium.precursorRange0 = new Range(Mass.MZFromMass(mTRAQMedium.experimentalMass, charge), new Tolerance(Form1.SEARCHPPM, ToleranceType.PPM));
                mTRAQMedium.precursorRange1 = new Range(Mass.MZFromMass(mTRAQMedium.experimentalMass + Constants.Neutron, charge), new Tolerance(Form1.SEARCHPPM, ToleranceType.PPM));
                mTRAQMedium.precursorRange2 = new Range(Mass.MZFromMass(mTRAQMedium.experimentalMass + (2.0 * Constants.Neutron), charge), new Tolerance(Form1.SEARCHPPM, ToleranceType.PPM));
                mTRAQMedium.precursorRange3 = new Range(Mass.MZFromMass(mTRAQMedium.experimentalMass + (3.0 * Constants.Neutron), charge), new Tolerance(Form1.SEARCHPPM, ToleranceType.PPM));

                mTRAQHeavy.calculatedIntraMZDifference = ((double)lysineCount * lysineVersion) / (double)charge;
                mTRAQHeavy.precursorRange0 = new Range(Mass.MZFromMass(mTRAQHeavy.experimentalMass, charge), new Tolerance(Form1.SEARCHPPM, ToleranceType.PPM));
                mTRAQHeavy.precursorRange1 = new Range(Mass.MZFromMass(mTRAQHeavy.experimentalMass + Constants.Neutron, charge), new Tolerance(Form1.SEARCHPPM, ToleranceType.PPM));
                mTRAQHeavy.precursorRange2 = new Range(Mass.MZFromMass(mTRAQHeavy.experimentalMass + (2.0 * Constants.Neutron), charge), new Tolerance(Form1.SEARCHPPM, ToleranceType.PPM));
                mTRAQHeavy.precursorRange3 = new Range(Mass.MZFromMass(mTRAQHeavy.experimentalMass + (3.0 * Constants.Neutron), charge), new Tolerance(Form1.SEARCHPPM, ToleranceType.PPM));               
            }
        }

        //Returns the peak with the maximum intensity within the specified m/z range
        public LabeledPeak findMaxPeak(Spectrum fullScanSpectrum, Range range)
        {
            Spectrum spectrum = new Spectrum();
            spectrum = fullScanSpectrum;
            LabeledPeak maxPeak = (LabeledPeak) spectrum.GetLargestPeak(range);
            /*
            bool rightCharge = false;
            while (!rightCharge && maxPeak != null)
            {
                if (maxPeak.Charge != 0 && maxPeak.Charge == chargeState)
                {
                    rightCharge = true;
                }
                else
                {
                    spectrum.Remove(new MZRange(maxPeak.MZ, new Tolerance(5.0, ToleranceType.PPM)));
                    maxPeak = (LabeledPeak) spectrum.GetLargestPeak(range);
                    Console.WriteLine("wrong charge for max peak");
                }
            }
            */

            if (maxPeak != null && maxPeak.SN > 2)
            {
                return maxPeak;
            }
            else
            {
                return null;
            }
        }
        
        /*Looks for HILAC pairs in monoisotope and first 3 isotopes
         *Once a HILAC pair is found in any isotopic version, a HILAC pair is created
         *Once a HILAC pair is created, light & heavy intensities are stored for the isotopic version that caused the creation and any subsequent isotopic versions
         * */
        public void findHILACPair(LabeledPeak maxPeak0, LabeledPeak maxPeak1, LabeledPeak maxPeak2, LabeledPeak maxPeak3, ThermoRawFileScan currentFullScan)
        {
            HILACPair mono = null;
            HILACPair first = null;
            HILACPair second = null;
            HILACPair third = null;
            
            mono = findHILACPair(charge, maxPeak0, currentFullScan);
            if (mono == null) //No HILAC pair found in monoisotope -- check 1st
            {
                first = findHILACPair(charge, maxPeak1, currentFullScan);
                if (first == null) //No HLAC pair found in first isotope -- check 2nd
                {
                    second= findHILACPair(charge, maxPeak2, currentFullScan);
                    if (second == null) //No HILAC pair found in first 2 isotopes -- check 3rd
                    {
                        third = findHILACPair(charge, maxPeak3, currentFullScan);
                        if (third == null) //No HILAC pair found in all 3 isotopes -- done
                        {
                        }
                        else //HILAC pair found in 3rd isotope -- append 3rd
                        {
                            appendHILACPair(charge, maxPeak3, currentFullScan, third, 3);
                        }
                    }
                    else //HILAC pair found in 2nd isotope -- append 2nd & 3rd
                    {
                        appendHILACPair(charge, maxPeak2, currentFullScan, second, 2);
                        appendHILACPair(charge, maxPeak3, currentFullScan, second, 3);
                    }
                }
                else //HILAC pair found in first isotope -- append 1st, 2nd & 3rd
                {
                    appendHILACPair(charge, maxPeak1, currentFullScan, first, 1);
                    appendHILACPair(charge, maxPeak2, currentFullScan, first, 2);
                    appendHILACPair(charge, maxPeak3, currentFullScan, first, 3);
                }
            }
            else //HILAC pair found in monoisotope -- append mono, 1st, 2nd & 3rd
            {
                appendHILACPair(charge, maxPeak0, currentFullScan, mono, 0);
                appendHILACPair(charge, maxPeak1, currentFullScan, mono, 1);
                appendHILACPair(charge, maxPeak2, currentFullScan, mono, 2);
                appendHILACPair(charge, maxPeak3, currentFullScan, mono, 3);
            }
        }
        
        //Returns a HILACPair associated with the specified peak, if one exists
        public HILACPair findHILACPair(int chargeState, LabeledPeak maxPeak, ThermoRawFileScan currentFullScan)
        {
            Spectrum fullScanSpectrum = currentFullScan.Spectrum;

            if (maxPeak != null)
            {
                LabeledPeak light = checkPartner(maxPeak, fullScanSpectrum, -1);
                LabeledPeak heavy = checkPartner(maxPeak, fullScanSpectrum, 1);
                HILACPair pair;

                if (light != null || heavy != null) //Either search return a peak within the tolerance
                {
                    pair = new HILACPair(this, currentFullScan.ScanNum);
                    hILACPairs.Add(pair);
                    return pair;
                }
                else //Neither search returns a peak within the tolerance & above intensity threshold
                {
                    return null;
                }
            }
            else //No "identification" peak found within range
            {
                return null;
            }
        }

        //Appends isotope information to an existing HILAC pair
        public void appendHILACPair(int chargeState, LabeledPeak maxPeak, ThermoRawFileScan currentFullScan, HILACPair pair, int isotope)
        {
            Spectrum fullScanSpectrum = currentFullScan.Spectrum;

            if (maxPeak != null)
            {
                LabeledPeak light = checkPartner(maxPeak, fullScanSpectrum, -1);
                LabeledPeak heavy = checkPartner(maxPeak, fullScanSpectrum, 1);

                if (light != null && heavy != null) //Both searches return a peak within the tolerance
                {
                    //Take most intense peak as "true" partner
                    if (light.Intensity >= heavy.Intensity) //"Light" is partner 
                    {
                        pair.isotopes[isotope] = new HILACPairIsotope(pair, isotope, light.Intensity, maxPeak.Intensity);
                    }

                    else //"Heavy" is partner
                    {
                        pair.isotopes[isotope] = new HILACPairIsotope(pair, isotope, maxPeak.Intensity, heavy.Intensity);
                    }
                }

                else if (light != null) //"Light" is partner ("heavy" is null)
                {
                    pair.isotopes[isotope] = new HILACPairIsotope(pair, isotope, light.Intensity, maxPeak.Intensity);
                }

                else if (heavy != null) //"Heavy" is partner ("light" is null)
                {
                    pair.isotopes[isotope] = new HILACPairIsotope(pair, isotope, maxPeak.Intensity, heavy.Intensity);
                }

                else //Neither search returns a peak within the tolerance & above intensity threshold
                {
                    if (reQuantify && isotope == 0 && pair.noiseBandCapped)
                    {
                        pair.isotopes[isotope] = new HILACPairIsotope(pair, isotope, maxPeak.Noise, maxPeak.Intensity);
                        pair.isotopes[isotope].noiseBandCapped = true;
                    }
                    else
                    {
                        pair.isotopes[isotope] = null;
                    }
                }
            }
            else //Max peak is null
            {
                pair.isotopes[isotope] = null;
            }
        }

        //Checks above and below specified peak for a partner peak -- returns null if a peak is not found above the noise level
        public LabeledPeak checkPartner(LabeledPeak maxPeak, Spectrum fullScanSpectrum, int partner)
        {
            Range potentialPartnerRange;
            LabeledPeak potentialPartnerPeak;
            if (partner < 0) //Look for light partner
            {
                potentialPartnerRange = new Range(maxPeak.MZ - calculatedIntraMZDifference, Form1.DELTAMZALLOWED);
            }
            else //Look for heavy partner
            {
                potentialPartnerRange = new Range(maxPeak.MZ + calculatedIntraMZDifference, Form1.DELTAMZALLOWED);         
            }
            potentialPartnerPeak = (LabeledPeak)fullScanSpectrum.GetLargestPeak(potentialPartnerRange);
            if (potentialPartnerPeak != null && potentialPartnerPeak.SN > 2)
            {
                return potentialPartnerPeak;
            }
            else
            {
                return null;
            }
        }

        //Calculate the median ratio for the peptide using a sorted list of quantifiable HILAC pairs
        public void calculateListMedian()
        {
            if ((double)(quantifiableHILACPairIsotopes.Count) % 2 != 0) //Odd # of HILAC pairs
            {
                int singleIndex = quantifiableHILACPairIsotopes.Count / 2;
                medianHeavyToLightRatio = quantifiableHILACPairIsotopes.Values[singleIndex].ratio;
            }
            else //Even # of HILAC Pairs
            {
                int doubleIndex1 = (quantifiableHILACPairIsotopes.Count / 2) - 1;
                int doubleIndex2 = doubleIndex1 + 1;
                medianHeavyToLightRatio = (quantifiableHILACPairIsotopes.Values[doubleIndex1].ratio + quantifiableHILACPairIsotopes.Values[doubleIndex2].ratio) / 2.0;
            }
        }

        //Quantify peptide -- null AveragePeptideID returned if quantifiable
        public bool quantify()
        {
            //No HILAC pairs
            if (hILACPairs.Count == 0)
            {
                reasonNQ = 1;
                if (Form1.ULTRAHILAC)
                {
                    parent.reasonNQ = 1;
                }
                return false;
            }
            //Fewer than 3 total HILAC pairs -- do not quantify
            else if (hILACPairs.Count < 3)
            {
                reasonNQ = 2;
                if (Form1.ULTRAHILAC)
                {
                    parent.reasonNQ = 2;
                }
                return false;
            }
            //More than 3 total HILAC pairs -- continue quantifying
            else
            {                
                //Align MS1 scans of HILAC pair components so that peak apexes are as close as possible (correction for retention time shifts)
                alignHILACPairs();

                //Filter out MS1 scans that are below quantitation level for either light or heavy
                if (alignedhILACPairs != null && alignedhILACPairs.Count > 2)
                {
                    filterHILACPairs();

                    //Assemble peptide quantitation
                    if (filteredhILACPairs != null && filteredhILACPairs.Count > 2)
                    {
                        quantifyHILACPairs();

                        if (quantifiableHILACPairIsotopes != null)
                        {
                            calculateListMedian();
                            return true;
                        }
                        else
                        {
                            reasonNQ = 5;
                            if (Form1.ULTRAHILAC)
                            {
                                parent.reasonNQ = 5;
                            }
                            return false;
                        }
                    }
                    else
                    {
                        reasonNQ = 4;
                        if (Form1.ULTRAHILAC)
                        {
                            parent.reasonNQ = 4;
                        }
                        return false;
                    }
                }
                else
                {
                    reasonNQ = 3;
                    if (Form1.ULTRAHILAC)
                    {
                        parent.reasonNQ = 3;
                    }
                    return false;
                }
            }
        }

        //Correct for retention time shifts between light and heavy by aligning peak apexes
        public void alignHILACPairs()
        {
            if (maxLight == null || maxHeavy == null)
            {
                alignedhILACPairs = null;
            }
            else
            {
                if (peakMaxShift == 0) //No alignment needed
                {
                    alignedhILACPairs = hILACPairs;
                }
                else if (peakMaxShift > 0) //Heavy elutes after light --> shift heavy down # scans
                {
                    alignedhILACPairs = hILACPairs;

                    //Align HILACPairs
                    for (int i = 0; i < alignedhILACPairs.Count - 1 - peakMaxShift; i++)
                    {
                        for (int j = 0; j < 4; j++)
                        {
                            if (hILACPairs.ElementAt(i).isotopes[j] != null && hILACPairs.ElementAt(i + peakMaxShift).isotopes[j] != null)
                            {
                                alignedhILACPairs.ElementAt(i).isotopes[j].heavy = hILACPairs.ElementAt(i + peakMaxShift).isotopes[j].heavy;
                            }
                        }
                    }

                    //Clean up light extra by removing last elements
                    for (int i = 0; i < peakMaxShift; i++)
                    {
                        alignedhILACPairs.RemoveAt(alignedhILACPairs.Count - 1);
                    }
                }
                else //Light elutes after heavy --> shift light down # scans
                {
                    int newPeakMaxShift = Math.Abs(peakMaxShift);
                    alignedhILACPairs = hILACPairs;

                    //Align HILACPairs
                    for (int i = 0; i < alignedhILACPairs.Count - 1 - newPeakMaxShift; i++)
                    {
                        for (int j = 0; j < 4; j++)
                        {
                            if (hILACPairs.ElementAt(i).isotopes[j] != null && hILACPairs.ElementAt(i + newPeakMaxShift).isotopes[j] != null)
                            {
                                alignedhILACPairs.ElementAt(i).isotopes[j].light = hILACPairs.ElementAt(i + newPeakMaxShift).isotopes[j].light;
                            }
                        }
                    }

                    //Clean up heavy extra by removing last elements
                    for (int i = 0; i < newPeakMaxShift; i++)
                    {
                        alignedhILACPairs.RemoveAt(alignedhILACPairs.Count - 1);
                    }
                }
            }
        }

        //Only include HILAC pairs for quantitation that have an isotope in which both light and heavy are above intensity cut off (1/e of the maximum)
        public void filterHILACPairs()
        {
            if (alignedhILACPairs != null)
            {
                filteredhILACPairs = new List<HILACPair>();
                bool goodIsotope;
                int counter;

                foreach (HILACPair pair in alignedhILACPairs)
                {
                    goodIsotope = false;
                    counter = 0;

                    while (!goodIsotope && counter < pair.isotopes.Count())
                    {
                        if (pair.isotopes[counter] != null && pair.isotopes[counter].light != null && pair.isotopes[counter].heavy != null)
                        {
                            if (pair.isotopes[counter].light.intensity >= lightIntensityCutOff || pair.isotopes[counter].heavy.intensity >= heavyIntensityCutOff)
                            {
                                goodIsotope = true;
                            }
                        }
                        counter++;
                    }

                    if (goodIsotope)
                    {
                        filteredhILACPairs.Add(pair);
                    }
                }
            }
        }

        public void quantifyHILACPairs()
        {
            quantifiableHILACPairIsotopes = new SortedList<double, HILACPairIsotope>();

            foreach (HILACPair pair in filteredhILACPairs)
            {
                foreach (HILACPairIsotope isotope in pair.isotopes)
                {
                    if (isotope != null)
                    {
                        try
                        {
                            quantifiableHILACPairIsotopes.Add(isotope.ratio, isotope);
                        }
                        catch (ArgumentException)
                        {
                            try
                            {
                                quantifiableHILACPairIsotopes.Add(isotope.ratio + 0.00000000001, isotope);
                            }
                            catch (ArgumentException)
                            {
                                try
                                {
                                    quantifiableHILACPairIsotopes.Add(isotope.ratio + 0.0000000001, isotope);
                                }
                                catch (ArgumentException)
                                {
                                    try
                                    {
                                        quantifiableHILACPairIsotopes.Add(isotope.ratio + 0.000000001, isotope);
                                    }
                                    catch (ArgumentException)
                                    {
                                        quantifiableHILACPairIsotopes.Add(isotope.ratio + 0.00000001, isotope);
                                    }
                                }
                            }
                        }
                    }
                }
                lightTotalIntensity += pair.totalLightIntensity;
                heavyTotalIntensity += pair.totalHeavyIntensity;
            }
                    
            XYPoint[] points = new XYPoint[quantifiableHILACPairIsotopes.Count];
            for (int i = 0; i < points.Count(); i++)
            {
                points[i] = new XYPoint(quantifiableHILACPairIsotopes.ElementAt(i).Value.light.intensity, quantifiableHILACPairIsotopes.ElementAt(i).Value.heavy.intensity);
            }

            double slope = 0.0;
            double yIntercept = 0.0;
            LeastSquaresFitLinear(points, points.Count(), ref slope, ref yIntercept);

            slopeHeavyToLightRatio = slope;
        }

        public static void LeastSquaresFitLinear(XYPoint[] points, int numPoints, ref double M, ref double B)
        {
            //Gives best fit of data to line Y = MC + B  
            double x1, y1, xy, x2, J;
            int i;

            x1 = 0.0;
            y1 = 0.0;
            xy = 0.0;
            x2 = 0.0;

            for (i = 0; i < numPoints; i++)
            {
                x1 = x1 + points[i].X;
                y1 = y1 + points[i].Y;
                xy = xy + points[i].X * points[i].Y;
                x2 = x2 + points[i].X * points[i].X;
            }

            J = ((double)numPoints * x2) - (x1 * x1);
            if (J != 0.0)
            {
                M = (((double)numPoints * xy) - (x1 * y1)) / J;
                M = System.Math.Floor(1.0E3 * M + 0.5) / 1.0E3;
                B = ((y1 * x2) - (x1 * xy)) / J;
                B = System.Math.Floor(1.0E3 * B + 0.5) / 1.0E3;
            }
            else
            {
                M = 0;
                B = 0;
            }
        }  

        public string getReasonNQ()
        {
            if (reasonNQ == 0)
            {
                return "# lysines: 0";
            }
            else if (reasonNQ == 1)
            {
                return "# quantifiable pairs: 0";
            }
            else if (reasonNQ == 2)
            {
                return "not enough quantifiable pairs";
            }
            else if (reasonNQ == 3)
            {
                return "not able to align pairs";
            }
            else if (reasonNQ == 4)
            {
                return "not able to intensity-filter aligned pairs";
            }
            else
            {
                return "not able to quantify intensity-filtered pairs";
            }
        }

        /*
        public void countNBC()
        {
            int count = 0;
            List<HILACPair> hILACPairsNBCRemoved;
            foreach (HILACPair pair in hILACPairs)
            {
                if (pair.noiseBandCapped)
                {
                    count++;
                }
            }

            if (count > 0)
            {
                noiseBandCapped = true;
            }
            else
            {
                noiseBandCapped = false;
            }

            
            if (hILACPairs.Count - count > 5)
            {
                hILACPairsNBCRemoved = new List<HILACPair>();
                foreach (HILACPair pair in hILACPairs)
                {
                    if (!pair.noiseBandCapped)
                    {
                        hILACPairsNBCRemoved.Add(pair);
                    }
                }
                hILACPairs = hILACPairsNBCRemoved;
                noiseBandCapped = false;
            }
            else
            {
                noiseBandCapped = true;
            }
        }
        */

        public bool checkMaxPeak(LabeledPeak maxPeak, int chargeState)
        {
            if (maxPeak == null)
            {
                return false;
            }
            
            bool isLight = false;
            double lightMass = theoreticalMass - (0.018 * lysineCount);
            double heavyMass = theoreticalMass + (0.018 * lysineCount);
            double checkLightError = Tolerance.GetError(Mass.MassFromMz(maxPeak.MZ, chargeState), lightMass, ToleranceType.PPM);
            double checkHeavyError = Tolerance.GetError(Mass.MassFromMz(maxPeak.MZ, chargeState), heavyMass, ToleranceType.PPM);

            if (checkLightError < checkHeavyError)
            {
                isLight = true;
            }
            return isLight;
        }

        public void reQuantifyHILACPair(LabeledPeak maxPeak0, LabeledPeak maxPeak1, LabeledPeak maxPeak2, LabeledPeak maxPeak3, ThermoRawFileScan currentFullScan)
        {
            bool maxPeak0Light = checkMaxPeak(maxPeak0, charge);
            bool maxPeak1Light = checkMaxPeak(maxPeak1, charge);
            bool maxPeak2Light = checkMaxPeak(maxPeak2, charge);
            bool maxPeak3Light = checkMaxPeak(maxPeak3, charge);

            HILACPair mono = null;
            HILACPair first = null;
            HILACPair second = null;
            HILACPair third = null;

            mono = reQuantifyHILACPair(charge, maxPeak0, currentFullScan, maxPeak0Light);
            if (mono == null) //No HILAC pair found in monoisotope -- check 1st
            {
                first = reQuantifyHILACPair(charge, maxPeak1, currentFullScan, maxPeak1Light);
                if (first == null) //No HLAC pair found in first isotope -- check 2nd
                {
                    second = reQuantifyHILACPair(charge, maxPeak2, currentFullScan, maxPeak2Light);
                    if (second == null) //No HILAC pair found in first 2 isotopes -- check 3rd
                    {
                        third = reQuantifyHILACPair(charge, maxPeak3, currentFullScan, maxPeak3Light);
                        if (third == null) //No HILAC pair found in all 3 isotopes -- done
                        {
                            if (reQuantify)
                            {
                                HILACPair pair = new HILACPair(this, currentFullScan.ScanNum);
                                hILACPairs.Add(pair);
                                pair.noiseBandCapped = true;
                                appendReQuantifiedHILACPair(charge, maxPeak0, currentFullScan, pair, maxPeak0Light, 0);
                            }
                        }
                        else //HILAC pair found in 3rd isotope -- append 3rd
                        {
                            appendHILACPair(charge, maxPeak3, currentFullScan, third, 3);
                        }
                    }
                    else //HILAC pair found in 2nd isotope -- append 2nd & 3rd
                    {
                        appendHILACPair(charge, maxPeak2, currentFullScan, second, 2);
                        appendHILACPair(charge, maxPeak3, currentFullScan, second, 3);
                    }
                }
                else //HILAC pair found in first isotope -- append 1st, 2nd & 3rd
                {
                    appendHILACPair(charge, maxPeak1, currentFullScan, first, 1);
                    appendHILACPair(charge, maxPeak2, currentFullScan, first, 2);
                    appendHILACPair(charge, maxPeak3, currentFullScan, first, 3);
                }
            }
            else //HILAC pair found in monoisotope -- append mono, 1st, 2nd & 3rd
            {
                appendHILACPair(charge, maxPeak0, currentFullScan, mono, 0);
                appendHILACPair(charge, maxPeak1, currentFullScan, mono, 1);
                appendHILACPair(charge, maxPeak2, currentFullScan, mono, 2);
                appendHILACPair(charge, maxPeak3, currentFullScan, mono, 3);
            }

        }

        //Returns HILACPair if found; returns null if no such HILACPair is found
        public HILACPair reQuantifyHILACPair(int chargeState, LabeledPeak maxPeak, ThermoRawFileScan currentFullScan, bool light)
        {
            Spectrum fullScanSpectrum = currentFullScan.Spectrum;

            if (maxPeak != null)
            {
                LabeledPeak partner;
                HILACPair pair;

                if (light) //Find heavy partner
                {
                    partner = checkPartner(maxPeak, fullScanSpectrum, 1);

                    if (partner != null) //Partner peak
                    {
                        pair = new HILACPair(this, currentFullScan.ScanNum);
                        hILACPairs.Add(pair);
                        return pair;
                    }
                    else //Noise peak
                    {
                        return null;
                    }

                }
                else //Find light partner
                {
                    partner = checkPartner(maxPeak, fullScanSpectrum, -1);

                    if (partner != null) //Partner peak
                    {
                        pair = new HILACPair(this, currentFullScan.ScanNum);
                        hILACPairs.Add(pair);
                        return pair;
                    }
                    else //Noise peak
                    {
                        return null;
                    }
                }
            }
            else //No peak found within range
            {
                return null;
            }
        }

        public void appendReQuantifiedHILACPair(int chargeState, LabeledPeak maxPeak, ThermoRawFileScan currentFullScan, HILACPair pair, bool light, int isotope)
        {
            Spectrum fullScanSpectrum = currentFullScan.Spectrum;

            if (maxPeak != null)
            {
                IPeak partner;

                if (light) //Find heavy partner
                {
                    partner = checkPartner(maxPeak, fullScanSpectrum, 1);
                    
                    if (partner != null && partner.Intensity >= maxPeak.Noise) //Partner peak
                    {
                        pair.isotopes[isotope] = new HILACPairIsotope(pair, isotope, maxPeak.Intensity, partner.Intensity);
                    }
                    else //Noise peak
                    {
                        if (isotope == 0)
                        {
                            pair.isotopes[isotope] = new HILACPairIsotope(pair, isotope, maxPeak.Intensity, maxPeak.Noise);
                            pair.isotopes[isotope].noiseBandCapped = true;
                        }
                        else
                        {
                            pair.isotopes[isotope] = null;
                        }
                    }
                }
                else //Find light partner
                {
                    partner = checkPartner(maxPeak, fullScanSpectrum, -1);

                    if (partner != null && partner.Intensity >= maxPeak.Noise) //Partner peak
                    {
                        pair.isotopes[isotope] = new HILACPairIsotope(pair, isotope, partner.Intensity, maxPeak.Intensity);
                    }
                    else //Noise peak
                    {
                        if (isotope == 0)
                        {
                            pair.isotopes[isotope] = new HILACPairIsotope(pair, isotope, partner.Intensity, maxPeak.Intensity);
                            pair.isotopes[isotope].noiseBandCapped = true;
                        }
                        else
                        {
                            pair.isotopes[isotope] = null;
                        }
                    }
                }
            }
            else //No peak found within range
            {
                
            }
        }
    }
}
