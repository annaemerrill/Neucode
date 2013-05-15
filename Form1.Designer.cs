namespace Coon.NeuQuant
{
    partial class Form1
    {
        /// <summary>
        /// Required designer variable.
        /// </summary>
        private System.ComponentModel.IContainer components = null;

        /// <summary>
        /// Clean up any resources being used.
        /// </summary>
        /// <param name="disposing">true if managed resources should be disposed; otherwise, false.</param>
        protected override void Dispose(bool disposing)
        {
            if (disposing && (components != null))
            {
                components.Dispose();
            }
            base.Dispose(disposing);
        }

        #region Windows Form Designer generated code

        /// <summary>
        /// Required method for Designer support - do not modify
        /// the contents of this method with the code editor.
        /// </summary>
        private void InitializeComponent()
        {
            this.label1 = new System.Windows.Forms.Label();
            this.csvInputBox = new System.Windows.Forms.TextBox();
            this.label2 = new System.Windows.Forms.Label();
            this.outputFolderBox = new System.Windows.Forms.TextBox();
            this.label3 = new System.Windows.Forms.Label();
            this.rtWindow = new System.Windows.Forms.NumericUpDown();
            this.label4 = new System.Windows.Forms.Label();
            this.RAWBrowse = new System.Windows.Forms.Button();
            this.CSVBrowse = new System.Windows.Forms.Button();
            this.OutputBrowse = new System.Windows.Forms.Button();
            this.start = new System.Windows.Forms.Button();
            this.browseOutputLocation = new System.Windows.Forms.FolderBrowserDialog();
            this.browseTargetInput = new System.Windows.Forms.OpenFileDialog();
            this.noiseBandCap = new System.Windows.Forms.CheckBox();
            this.coalescence = new System.Windows.Forms.CheckBox();
            this.browseRawLocation = new System.Windows.Forms.FolderBrowserDialog();
            this.rawFileBox = new System.Windows.Forms.TextBox();
            this.signalToNoiseThreshold = new System.Windows.Forms.NumericUpDown();
            this.label6 = new System.Windows.Forms.Label();
            this.label5 = new System.Windows.Forms.Label();
            this.Isotopes = new System.Windows.Forms.NumericUpDown();
            this.Conversion = new System.Windows.Forms.CheckBox();
            this.label7 = new System.Windows.Forms.Label();
            this.filterProfiles = new System.Windows.Forms.CheckBox();
            this.label8 = new System.Windows.Forms.Label();
            this.label9 = new System.Windows.Forms.Label();
            this.SILACLys8CN = new System.Windows.Forms.RadioButton();
            this.SILACLys8D = new System.Windows.Forms.RadioButton();
            this.NeuCodeLys1 = new System.Windows.Forms.RadioButton();
            this.NeuCodeLys8Duplex = new System.Windows.Forms.RadioButton();
            this.NeuCodeLys8Triplex = new System.Windows.Forms.RadioButton();
            this.NeuCodeLys8Fourplex = new System.Windows.Forms.RadioButton();
            this.NeuCodeLys8Sixplex = new System.Windows.Forms.RadioButton();
            this.label10 = new System.Windows.Forms.Label();
            this.NeuCodeLeu7Duplex = new System.Windows.Forms.RadioButton();
            this.SILACLeu7CN = new System.Windows.Forms.RadioButton();
            this.SILACLeu7D = new System.Windows.Forms.RadioButton();
            this.CarbamylCN = new System.Windows.Forms.RadioButton();
            this.FourplexL = new System.Windows.Forms.RadioButton();
            this.FourplexM = new System.Windows.Forms.RadioButton();
            this.FourplexH = new System.Windows.Forms.RadioButton();
            this.Twelveplex = new System.Windows.Forms.RadioButton();
            this.mTRAQ = new System.Windows.Forms.RadioButton();
            this.Arg = new System.Windows.Forms.RadioButton();
            this.Leu = new System.Windows.Forms.RadioButton();
            this.IncompleteIncorporation = new System.Windows.Forms.CheckBox();
            this.label11 = new System.Windows.Forms.Label();
            this.label12 = new System.Windows.Forms.Label();
            this.PeakSeparation = new System.Windows.Forms.NumericUpDown();
            this.QuantResolution = new System.Windows.Forms.NumericUpDown();
            ((System.ComponentModel.ISupportInitialize)(this.rtWindow)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.signalToNoiseThreshold)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.Isotopes)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.PeakSeparation)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.QuantResolution)).BeginInit();
            this.SuspendLayout();
            // 
            // label1
            // 
            this.label1.AutoSize = true;
            this.label1.Location = new System.Drawing.Point(12, 30);
            this.label1.Name = "label1";
            this.label1.Size = new System.Drawing.Size(116, 13);
            this.label1.TabIndex = 1;
            this.label1.Text = "Location of .RAW Files";
            // 
            // csvInputBox
            // 
            this.csvInputBox.Location = new System.Drawing.Point(15, 90);
            this.csvInputBox.Name = "csvInputBox";
            this.csvInputBox.Size = new System.Drawing.Size(449, 20);
            this.csvInputBox.TabIndex = 2;
            // 
            // label2
            // 
            this.label2.AutoSize = true;
            this.label2.Location = new System.Drawing.Point(9, 74);
            this.label2.Name = "label2";
            this.label2.Size = new System.Drawing.Size(73, 13);
            this.label2.TabIndex = 3;
            this.label2.Text = "target.csv File";
            // 
            // outputFolderBox
            // 
            this.outputFolderBox.Location = new System.Drawing.Point(12, 145);
            this.outputFolderBox.Name = "outputFolderBox";
            this.outputFolderBox.Size = new System.Drawing.Size(452, 20);
            this.outputFolderBox.TabIndex = 4;
            // 
            // label3
            // 
            this.label3.AutoSize = true;
            this.label3.Location = new System.Drawing.Point(12, 129);
            this.label3.Name = "label3";
            this.label3.Size = new System.Drawing.Size(102, 13);
            this.label3.TabIndex = 5;
            this.label3.Text = "Output File Location";
            // 
            // rtWindow
            // 
            this.rtWindow.DecimalPlaces = 2;
            this.rtWindow.Increment = new decimal(new int[] {
            25,
            0,
            0,
            131072});
            this.rtWindow.Location = new System.Drawing.Point(192, 338);
            this.rtWindow.Maximum = new decimal(new int[] {
            5,
            0,
            0,
            0});
            this.rtWindow.Name = "rtWindow";
            this.rtWindow.Size = new System.Drawing.Size(47, 20);
            this.rtWindow.TabIndex = 6;
            this.rtWindow.Value = new decimal(new int[] {
            50,
            0,
            0,
            131072});
            // 
            // label4
            // 
            this.label4.AutoSize = true;
            this.label4.Location = new System.Drawing.Point(189, 322);
            this.label4.Name = "label4";
            this.label4.Size = new System.Drawing.Size(47, 13);
            this.label4.TabIndex = 7;
            this.label4.Text = "RT (min)";
            // 
            // RAWBrowse
            // 
            this.RAWBrowse.Location = new System.Drawing.Point(491, 46);
            this.RAWBrowse.Name = "RAWBrowse";
            this.RAWBrowse.Size = new System.Drawing.Size(75, 23);
            this.RAWBrowse.TabIndex = 12;
            this.RAWBrowse.Text = "Browse";
            this.RAWBrowse.UseVisualStyleBackColor = true;
            this.RAWBrowse.Click += new System.EventHandler(this.RAWBrowse_Click);
            // 
            // CSVBrowse
            // 
            this.CSVBrowse.Location = new System.Drawing.Point(491, 90);
            this.CSVBrowse.Name = "CSVBrowse";
            this.CSVBrowse.Size = new System.Drawing.Size(75, 23);
            this.CSVBrowse.TabIndex = 13;
            this.CSVBrowse.Text = "Browse";
            this.CSVBrowse.UseVisualStyleBackColor = true;
            this.CSVBrowse.Click += new System.EventHandler(this.CSVBrowse_Click);
            // 
            // OutputBrowse
            // 
            this.OutputBrowse.Location = new System.Drawing.Point(491, 143);
            this.OutputBrowse.Name = "OutputBrowse";
            this.OutputBrowse.Size = new System.Drawing.Size(75, 23);
            this.OutputBrowse.TabIndex = 14;
            this.OutputBrowse.Text = "Browse";
            this.OutputBrowse.UseVisualStyleBackColor = true;
            this.OutputBrowse.Click += new System.EventHandler(this.OutputBrowse_Click);
            // 
            // start
            // 
            this.start.Location = new System.Drawing.Point(491, 208);
            this.start.Name = "start";
            this.start.Size = new System.Drawing.Size(75, 23);
            this.start.TabIndex = 17;
            this.start.Text = "Go!";
            this.start.UseVisualStyleBackColor = true;
            this.start.Click += new System.EventHandler(this.start_Click);
            // 
            // browseTargetInput
            // 
            this.browseTargetInput.FileName = "openFileDialog1";
            // 
            // noiseBandCap
            // 
            this.noiseBandCap.AutoSize = true;
            this.noiseBandCap.Location = new System.Drawing.Point(393, 367);
            this.noiseBandCap.Name = "noiseBandCap";
            this.noiseBandCap.Size = new System.Drawing.Size(183, 17);
            this.noiseBandCap.TabIndex = 18;
            this.noiseBandCap.Text = "Noise Band Cap Missing Channel";
            this.noiseBandCap.UseVisualStyleBackColor = true;
            // 
            // coalescence
            // 
            this.coalescence.AutoSize = true;
            this.coalescence.Location = new System.Drawing.Point(393, 390);
            this.coalescence.Name = "coalescence";
            this.coalescence.Size = new System.Drawing.Size(147, 17);
            this.coalescence.TabIndex = 19;
            this.coalescence.Text = "Track Peak Coalescence";
            this.coalescence.UseVisualStyleBackColor = true;
            // 
            // rawFileBox
            // 
            this.rawFileBox.Location = new System.Drawing.Point(15, 49);
            this.rawFileBox.Name = "rawFileBox";
            this.rawFileBox.Size = new System.Drawing.Size(449, 20);
            this.rawFileBox.TabIndex = 20;
            // 
            // signalToNoiseThreshold
            // 
            this.signalToNoiseThreshold.DecimalPlaces = 1;
            this.signalToNoiseThreshold.Increment = new decimal(new int[] {
            5,
            0,
            0,
            65536});
            this.signalToNoiseThreshold.Location = new System.Drawing.Point(255, 338);
            this.signalToNoiseThreshold.Maximum = new decimal(new int[] {
            10,
            0,
            0,
            0});
            this.signalToNoiseThreshold.Name = "signalToNoiseThreshold";
            this.signalToNoiseThreshold.Size = new System.Drawing.Size(38, 20);
            this.signalToNoiseThreshold.TabIndex = 21;
            this.signalToNoiseThreshold.Value = new decimal(new int[] {
            30,
            0,
            0,
            65536});
            // 
            // label6
            // 
            this.label6.AutoSize = true;
            this.label6.Location = new System.Drawing.Point(252, 322);
            this.label6.Name = "label6";
            this.label6.Size = new System.Drawing.Size(50, 13);
            this.label6.TabIndex = 22;
            this.label6.Text = "Min. S/N";
            this.label6.Click += new System.EventHandler(this.label6_Click);
            // 
            // label5
            // 
            this.label5.AutoSize = true;
            this.label5.Location = new System.Drawing.Point(324, 324);
            this.label5.Name = "label5";
            this.label5.Size = new System.Drawing.Size(47, 13);
            this.label5.TabIndex = 23;
            this.label5.Text = "Isotopes";
            // 
            // Isotopes
            // 
            this.Isotopes.Location = new System.Drawing.Point(327, 340);
            this.Isotopes.Maximum = new decimal(new int[] {
            10,
            0,
            0,
            0});
            this.Isotopes.Minimum = new decimal(new int[] {
            1,
            0,
            0,
            0});
            this.Isotopes.Name = "Isotopes";
            this.Isotopes.Size = new System.Drawing.Size(30, 20);
            this.Isotopes.TabIndex = 24;
            this.Isotopes.Value = new decimal(new int[] {
            3,
            0,
            0,
            0});
            // 
            // Conversion
            // 
            this.Conversion.AutoSize = true;
            this.Conversion.Location = new System.Drawing.Point(393, 344);
            this.Conversion.Name = "Conversion";
            this.Conversion.Size = new System.Drawing.Size(157, 17);
            this.Conversion.TabIndex = 25;
            this.Conversion.Text = "Check for Label Conversion";
            this.Conversion.UseVisualStyleBackColor = true;
            // 
            // label7
            // 
            this.label7.AutoSize = true;
            this.label7.Location = new System.Drawing.Point(390, 296);
            this.label7.Name = "label7";
            this.label7.Size = new System.Drawing.Size(111, 13);
            this.label7.TabIndex = 26;
            this.label7.Text = "Extra Analysis Options";
            // 
            // filterProfiles
            // 
            this.filterProfiles.AutoSize = true;
            this.filterProfiles.Location = new System.Drawing.Point(393, 321);
            this.filterProfiles.Name = "filterProfiles";
            this.filterProfiles.Size = new System.Drawing.Size(152, 17);
            this.filterProfiles.TabIndex = 27;
            this.filterProfiles.Text = "Filter Quant Elution Profiles";
            this.filterProfiles.UseVisualStyleBackColor = true;
            // 
            // label8
            // 
            this.label8.AutoSize = true;
            this.label8.Location = new System.Drawing.Point(12, 179);
            this.label8.Name = "label8";
            this.label8.Size = new System.Drawing.Size(96, 13);
            this.label8.TabIndex = 28;
            this.label8.Text = "Metabolic Labeling";
            // 
            // label9
            // 
            this.label9.AutoSize = true;
            this.label9.Location = new System.Drawing.Point(219, 179);
            this.label9.Name = "label9";
            this.label9.Size = new System.Drawing.Size(93, 13);
            this.label9.TabIndex = 29;
            this.label9.Text = "Chemical Labeling";
            // 
            // SILACLys8CN
            // 
            this.SILACLys8CN.AutoSize = true;
            this.SILACLys8CN.Location = new System.Drawing.Point(15, 208);
            this.SILACLys8CN.Name = "SILACLys8CN";
            this.SILACLys8CN.Size = new System.Drawing.Size(113, 17);
            this.SILACLys8CN.TabIndex = 30;
            this.SILACLys8CN.TabStop = true;
            this.SILACLys8CN.Text = "Lys SILAC: +8(CN)";
            this.SILACLys8CN.UseVisualStyleBackColor = true;
            // 
            // SILACLys8D
            // 
            this.SILACLys8D.AutoSize = true;
            this.SILACLys8D.Location = new System.Drawing.Point(15, 231);
            this.SILACLys8D.Name = "SILACLys8D";
            this.SILACLys8D.Size = new System.Drawing.Size(106, 17);
            this.SILACLys8D.TabIndex = 31;
            this.SILACLys8D.TabStop = true;
            this.SILACLys8D.Text = "Lys SILAC: +8(D)";
            this.SILACLys8D.UseVisualStyleBackColor = true;
            // 
            // NeuCodeLys1
            // 
            this.NeuCodeLys1.AutoSize = true;
            this.NeuCodeLys1.Location = new System.Drawing.Point(15, 297);
            this.NeuCodeLys1.Name = "NeuCodeLys1";
            this.NeuCodeLys1.Size = new System.Drawing.Size(141, 17);
            this.NeuCodeLys1.TabIndex = 32;
            this.NeuCodeLys1.TabStop = true;
            this.NeuCodeLys1.Text = "Lys NeuCode: +1(6mDa)";
            this.NeuCodeLys1.UseVisualStyleBackColor = true;
            // 
            // NeuCodeLys8Duplex
            // 
            this.NeuCodeLys8Duplex.AutoSize = true;
            this.NeuCodeLys8Duplex.Location = new System.Drawing.Point(15, 320);
            this.NeuCodeLys8Duplex.Name = "NeuCodeLys8Duplex";
            this.NeuCodeLys8Duplex.Size = new System.Drawing.Size(147, 17);
            this.NeuCodeLys8Duplex.TabIndex = 33;
            this.NeuCodeLys8Duplex.TabStop = true;
            this.NeuCodeLys8Duplex.Text = "Lys NeuCode: +8(36mDa)";
            this.NeuCodeLys8Duplex.UseVisualStyleBackColor = true;
            // 
            // NeuCodeLys8Triplex
            // 
            this.NeuCodeLys8Triplex.AutoSize = true;
            this.NeuCodeLys8Triplex.Location = new System.Drawing.Point(15, 343);
            this.NeuCodeLys8Triplex.Name = "NeuCodeLys8Triplex";
            this.NeuCodeLys8Triplex.Size = new System.Drawing.Size(147, 17);
            this.NeuCodeLys8Triplex.TabIndex = 34;
            this.NeuCodeLys8Triplex.TabStop = true;
            this.NeuCodeLys8Triplex.Text = "Lys NeuCode: +8(18mDa)";
            this.NeuCodeLys8Triplex.UseVisualStyleBackColor = true;
            // 
            // NeuCodeLys8Fourplex
            // 
            this.NeuCodeLys8Fourplex.AutoSize = true;
            this.NeuCodeLys8Fourplex.Location = new System.Drawing.Point(15, 366);
            this.NeuCodeLys8Fourplex.Name = "NeuCodeLys8Fourplex";
            this.NeuCodeLys8Fourplex.Size = new System.Drawing.Size(147, 17);
            this.NeuCodeLys8Fourplex.TabIndex = 35;
            this.NeuCodeLys8Fourplex.TabStop = true;
            this.NeuCodeLys8Fourplex.Text = "Lys NeuCode: +8(12mDa)";
            this.NeuCodeLys8Fourplex.UseVisualStyleBackColor = true;
            // 
            // NeuCodeLys8Sixplex
            // 
            this.NeuCodeLys8Sixplex.AutoSize = true;
            this.NeuCodeLys8Sixplex.Location = new System.Drawing.Point(15, 389);
            this.NeuCodeLys8Sixplex.Name = "NeuCodeLys8Sixplex";
            this.NeuCodeLys8Sixplex.Size = new System.Drawing.Size(141, 17);
            this.NeuCodeLys8Sixplex.TabIndex = 36;
            this.NeuCodeLys8Sixplex.TabStop = true;
            this.NeuCodeLys8Sixplex.Text = "Lys NeuCode: +8(6mDa)";
            this.NeuCodeLys8Sixplex.UseVisualStyleBackColor = true;
            // 
            // label10
            // 
            this.label10.AutoSize = true;
            this.label10.Location = new System.Drawing.Point(390, 179);
            this.label10.Name = "label10";
            this.label10.Size = new System.Drawing.Size(73, 13);
            this.label10.TabIndex = 37;
            this.label10.Text = "Cluster Labels";
            // 
            // NeuCodeLeu7Duplex
            // 
            this.NeuCodeLeu7Duplex.AutoSize = true;
            this.NeuCodeLeu7Duplex.Location = new System.Drawing.Point(15, 412);
            this.NeuCodeLeu7Duplex.Name = "NeuCodeLeu7Duplex";
            this.NeuCodeLeu7Duplex.Size = new System.Drawing.Size(149, 17);
            this.NeuCodeLeu7Duplex.TabIndex = 38;
            this.NeuCodeLeu7Duplex.TabStop = true;
            this.NeuCodeLeu7Duplex.Text = "Leu NeuCode: +7(18mDa)";
            this.NeuCodeLeu7Duplex.UseVisualStyleBackColor = true;
            // 
            // SILACLeu7CN
            // 
            this.SILACLeu7CN.AutoSize = true;
            this.SILACLeu7CN.Location = new System.Drawing.Point(15, 254);
            this.SILACLeu7CN.Name = "SILACLeu7CN";
            this.SILACLeu7CN.Size = new System.Drawing.Size(115, 17);
            this.SILACLeu7CN.TabIndex = 39;
            this.SILACLeu7CN.TabStop = true;
            this.SILACLeu7CN.Text = "Leu SILAC: +7(CN)";
            this.SILACLeu7CN.UseVisualStyleBackColor = true;
            // 
            // SILACLeu7D
            // 
            this.SILACLeu7D.AutoSize = true;
            this.SILACLeu7D.Location = new System.Drawing.Point(15, 277);
            this.SILACLeu7D.Name = "SILACLeu7D";
            this.SILACLeu7D.Size = new System.Drawing.Size(108, 17);
            this.SILACLeu7D.TabIndex = 40;
            this.SILACLeu7D.TabStop = true;
            this.SILACLeu7D.Text = "Leu SILAC: +7(D)";
            this.SILACLeu7D.UseVisualStyleBackColor = true;
            // 
            // CarbamylCN
            // 
            this.CarbamylCN.AutoSize = true;
            this.CarbamylCN.Location = new System.Drawing.Point(222, 204);
            this.CarbamylCN.Name = "CarbamylCN";
            this.CarbamylCN.Size = new System.Drawing.Size(68, 17);
            this.CarbamylCN.TabIndex = 41;
            this.CarbamylCN.TabStop = true;
            this.CarbamylCN.Text = "Carbamyl";
            this.CarbamylCN.UseVisualStyleBackColor = true;
            // 
            // FourplexL
            // 
            this.FourplexL.AutoSize = true;
            this.FourplexL.Location = new System.Drawing.Point(222, 227);
            this.FourplexL.Name = "FourplexL";
            this.FourplexL.Size = new System.Drawing.Size(79, 17);
            this.FourplexL.TabIndex = 42;
            this.FourplexL.TabStop = true;
            this.FourplexL.Text = "4plex: Light";
            this.FourplexL.UseVisualStyleBackColor = true;
            // 
            // FourplexM
            // 
            this.FourplexM.AutoSize = true;
            this.FourplexM.Location = new System.Drawing.Point(222, 250);
            this.FourplexM.Name = "FourplexM";
            this.FourplexM.Size = new System.Drawing.Size(93, 17);
            this.FourplexM.TabIndex = 43;
            this.FourplexM.TabStop = true;
            this.FourplexM.Text = "4plex: Medium";
            this.FourplexM.UseVisualStyleBackColor = true;
            // 
            // FourplexH
            // 
            this.FourplexH.AutoSize = true;
            this.FourplexH.Location = new System.Drawing.Point(222, 273);
            this.FourplexH.Name = "FourplexH";
            this.FourplexH.Size = new System.Drawing.Size(87, 17);
            this.FourplexH.TabIndex = 44;
            this.FourplexH.TabStop = true;
            this.FourplexH.Text = "4plex: Heavy";
            this.FourplexH.UseVisualStyleBackColor = true;
            // 
            // Twelveplex
            // 
            this.Twelveplex.AutoSize = true;
            this.Twelveplex.Location = new System.Drawing.Point(222, 296);
            this.Twelveplex.Name = "Twelveplex";
            this.Twelveplex.Size = new System.Drawing.Size(56, 17);
            this.Twelveplex.TabIndex = 45;
            this.Twelveplex.TabStop = true;
            this.Twelveplex.Text = "12plex";
            this.Twelveplex.UseVisualStyleBackColor = true;
            // 
            // mTRAQ
            // 
            this.mTRAQ.AutoSize = true;
            this.mTRAQ.Location = new System.Drawing.Point(393, 204);
            this.mTRAQ.Name = "mTRAQ";
            this.mTRAQ.Size = new System.Drawing.Size(63, 17);
            this.mTRAQ.TabIndex = 46;
            this.mTRAQ.TabStop = true;
            this.mTRAQ.Text = "mTRAQ";
            this.mTRAQ.UseVisualStyleBackColor = true;
            // 
            // Arg
            // 
            this.Arg.AutoSize = true;
            this.Arg.Location = new System.Drawing.Point(393, 227);
            this.Arg.Name = "Arg";
            this.Arg.Size = new System.Drawing.Size(41, 17);
            this.Arg.TabIndex = 47;
            this.Arg.TabStop = true;
            this.Arg.Text = "Arg";
            this.Arg.UseVisualStyleBackColor = true;
            // 
            // Leu
            // 
            this.Leu.AutoSize = true;
            this.Leu.Location = new System.Drawing.Point(393, 250);
            this.Leu.Name = "Leu";
            this.Leu.Size = new System.Drawing.Size(43, 17);
            this.Leu.TabIndex = 48;
            this.Leu.TabStop = true;
            this.Leu.Text = "Leu";
            this.Leu.UseVisualStyleBackColor = true;
            this.Leu.CheckedChanged += new System.EventHandler(this.NeuCodeLysLeu_CheckedChanged);
            // 
            // IncompleteIncorporation
            // 
            this.IncompleteIncorporation.AutoSize = true;
            this.IncompleteIncorporation.Location = new System.Drawing.Point(393, 413);
            this.IncompleteIncorporation.Name = "IncompleteIncorporation";
            this.IncompleteIncorporation.Size = new System.Drawing.Size(169, 17);
            this.IncompleteIncorporation.TabIndex = 49;
            this.IncompleteIncorporation.Text = "Check for Partial Incorporation";
            this.IncompleteIncorporation.UseVisualStyleBackColor = true;
            // 
            // label11
            // 
            this.label11.AutoSize = true;
            this.label11.Location = new System.Drawing.Point(201, 371);
            this.label11.Name = "label11";
            this.label11.Size = new System.Drawing.Size(38, 13);
            this.label11.TabIndex = 50;
            this.label11.Text = "FWxM";
            // 
            // label12
            // 
            this.label12.AutoSize = true;
            this.label12.Location = new System.Drawing.Point(284, 371);
            this.label12.Name = "label12";
            this.label12.Size = new System.Drawing.Size(73, 13);
            this.label12.TabIndex = 51;
            this.label12.Text = "Resolution (K)";
            // 
            // PeakSeparation
            // 
            this.PeakSeparation.Location = new System.Drawing.Point(204, 390);
            this.PeakSeparation.Minimum = new decimal(new int[] {
            1,
            0,
            0,
            0});
            this.PeakSeparation.Name = "PeakSeparation";
            this.PeakSeparation.Size = new System.Drawing.Size(35, 20);
            this.PeakSeparation.TabIndex = 52;
            this.PeakSeparation.Value = new decimal(new int[] {
            10,
            0,
            0,
            0});
            // 
            // QuantResolution
            // 
            this.QuantResolution.Increment = new decimal(new int[] {
            10,
            0,
            0,
            0});
            this.QuantResolution.Location = new System.Drawing.Point(287, 390);
            this.QuantResolution.Maximum = new decimal(new int[] {
            960,
            0,
            0,
            0});
            this.QuantResolution.Minimum = new decimal(new int[] {
            30,
            0,
            0,
            0});
            this.QuantResolution.Name = "QuantResolution";
            this.QuantResolution.Size = new System.Drawing.Size(44, 20);
            this.QuantResolution.TabIndex = 53;
            this.QuantResolution.Value = new decimal(new int[] {
            480,
            0,
            0,
            0});
            // 
            // Form1
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 13F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(589, 450);
            this.Controls.Add(this.QuantResolution);
            this.Controls.Add(this.PeakSeparation);
            this.Controls.Add(this.label12);
            this.Controls.Add(this.label11);
            this.Controls.Add(this.IncompleteIncorporation);
            this.Controls.Add(this.Leu);
            this.Controls.Add(this.Arg);
            this.Controls.Add(this.mTRAQ);
            this.Controls.Add(this.Twelveplex);
            this.Controls.Add(this.FourplexH);
            this.Controls.Add(this.FourplexM);
            this.Controls.Add(this.FourplexL);
            this.Controls.Add(this.CarbamylCN);
            this.Controls.Add(this.SILACLeu7D);
            this.Controls.Add(this.SILACLeu7CN);
            this.Controls.Add(this.NeuCodeLeu7Duplex);
            this.Controls.Add(this.label10);
            this.Controls.Add(this.NeuCodeLys8Sixplex);
            this.Controls.Add(this.NeuCodeLys8Fourplex);
            this.Controls.Add(this.NeuCodeLys8Triplex);
            this.Controls.Add(this.NeuCodeLys8Duplex);
            this.Controls.Add(this.NeuCodeLys1);
            this.Controls.Add(this.SILACLys8D);
            this.Controls.Add(this.SILACLys8CN);
            this.Controls.Add(this.label9);
            this.Controls.Add(this.label8);
            this.Controls.Add(this.filterProfiles);
            this.Controls.Add(this.label7);
            this.Controls.Add(this.Conversion);
            this.Controls.Add(this.Isotopes);
            this.Controls.Add(this.label5);
            this.Controls.Add(this.label6);
            this.Controls.Add(this.signalToNoiseThreshold);
            this.Controls.Add(this.rawFileBox);
            this.Controls.Add(this.coalescence);
            this.Controls.Add(this.noiseBandCap);
            this.Controls.Add(this.start);
            this.Controls.Add(this.OutputBrowse);
            this.Controls.Add(this.CSVBrowse);
            this.Controls.Add(this.RAWBrowse);
            this.Controls.Add(this.label4);
            this.Controls.Add(this.rtWindow);
            this.Controls.Add(this.label3);
            this.Controls.Add(this.outputFolderBox);
            this.Controls.Add(this.label2);
            this.Controls.Add(this.csvInputBox);
            this.Controls.Add(this.label1);
            this.Name = "Form1";
            this.Text = "Form1";
            this.Load += new System.EventHandler(this.Form1_Load);
            ((System.ComponentModel.ISupportInitialize)(this.rtWindow)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.signalToNoiseThreshold)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.Isotopes)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.PeakSeparation)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.QuantResolution)).EndInit();
            this.ResumeLayout(false);
            this.PerformLayout();

        }

        #endregion

        private System.Windows.Forms.Label label1;
        private System.Windows.Forms.TextBox csvInputBox;
        private System.Windows.Forms.Label label2;
        private System.Windows.Forms.TextBox outputFolderBox;
        private System.Windows.Forms.Label label3;
        private System.Windows.Forms.NumericUpDown rtWindow;
        private System.Windows.Forms.Label label4;
        private System.Windows.Forms.Button RAWBrowse;
        private System.Windows.Forms.Button CSVBrowse;
        private System.Windows.Forms.Button OutputBrowse;
        private System.Windows.Forms.Button start;
        private System.Windows.Forms.FolderBrowserDialog browseOutputLocation;
        private System.Windows.Forms.OpenFileDialog browseTargetInput;
        private System.Windows.Forms.CheckBox noiseBandCap;
        private System.Windows.Forms.CheckBox coalescence;
        private System.Windows.Forms.FolderBrowserDialog browseRawLocation;
        private System.Windows.Forms.TextBox rawFileBox;
        private System.Windows.Forms.NumericUpDown signalToNoiseThreshold;
        private System.Windows.Forms.Label label6;
        private System.Windows.Forms.Label label5;
        private System.Windows.Forms.NumericUpDown Isotopes;
        private System.Windows.Forms.CheckBox Conversion;
        private System.Windows.Forms.Label label7;
        private System.Windows.Forms.CheckBox filterProfiles;
        private System.Windows.Forms.Label label8;
        private System.Windows.Forms.Label label9;
        private System.Windows.Forms.RadioButton SILACLys8CN;
        private System.Windows.Forms.RadioButton SILACLys8D;
        private System.Windows.Forms.RadioButton NeuCodeLys1;
        private System.Windows.Forms.RadioButton NeuCodeLys8Duplex;
        private System.Windows.Forms.RadioButton NeuCodeLys8Triplex;
        private System.Windows.Forms.RadioButton NeuCodeLys8Fourplex;
        private System.Windows.Forms.RadioButton NeuCodeLys8Sixplex;
        private System.Windows.Forms.Label label10;
        private System.Windows.Forms.RadioButton NeuCodeLeu7Duplex;
        private System.Windows.Forms.RadioButton SILACLeu7CN;
        private System.Windows.Forms.RadioButton SILACLeu7D;
        private System.Windows.Forms.RadioButton CarbamylCN;
        private System.Windows.Forms.RadioButton FourplexL;
        private System.Windows.Forms.RadioButton FourplexM;
        private System.Windows.Forms.RadioButton FourplexH;
        private System.Windows.Forms.RadioButton Twelveplex;
        private System.Windows.Forms.RadioButton mTRAQ;
        private System.Windows.Forms.RadioButton Arg;
        private System.Windows.Forms.RadioButton Leu;
        private System.Windows.Forms.CheckBox IncompleteIncorporation;
        private System.Windows.Forms.Label label11;
        private System.Windows.Forms.Label label12;
        private System.Windows.Forms.NumericUpDown PeakSeparation;
        private System.Windows.Forms.NumericUpDown QuantResolution;
    }
}

