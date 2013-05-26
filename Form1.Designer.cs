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
            this.IncompleteIncorporation = new System.Windows.Forms.CheckBox();
            this.label11 = new System.Windows.Forms.Label();
            this.label12 = new System.Windows.Forms.Label();
            this.PeakSeparation = new System.Windows.Forms.NumericUpDown();
            this.QuantResolution = new System.Windows.Forms.NumericUpDown();
            this.mTRAQ = new System.Windows.Forms.CheckBox();
            this.Arg = new System.Windows.Forms.CheckBox();
            this.Leu = new System.Windows.Forms.CheckBox();
            this.label13 = new System.Windows.Forms.Label();
            this.MultipleInjections = new System.Windows.Forms.CheckBox();
            this.AGCBins = new System.Windows.Forms.CheckBox();
            this.Icat = new System.Windows.Forms.RadioButton();
            this.listBox1 = new System.Windows.Forms.ListBox();
            this.splitContainer1 = new System.Windows.Forms.SplitContainer();
            this.richTextBox1 = new System.Windows.Forms.RichTextBox();
            this.button1 = new System.Windows.Forms.Button();
            ((System.ComponentModel.ISupportInitialize)(this.rtWindow)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.signalToNoiseThreshold)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.Isotopes)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.PeakSeparation)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.QuantResolution)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.splitContainer1)).BeginInit();
            this.splitContainer1.Panel1.SuspendLayout();
            this.splitContainer1.Panel2.SuspendLayout();
            this.splitContainer1.SuspendLayout();
            this.SuspendLayout();
            // 
            // label1
            // 
            this.label1.AutoSize = true;
            this.label1.Location = new System.Drawing.Point(6, 77);
            this.label1.Name = "label1";
            this.label1.Size = new System.Drawing.Size(116, 13);
            this.label1.TabIndex = 1;
            this.label1.Text = "Location of .RAW Files";
            // 
            // label2
            // 
            this.label2.AutoSize = true;
            this.label2.Location = new System.Drawing.Point(6, 5);
            this.label2.Name = "label2";
            this.label2.Size = new System.Drawing.Size(73, 13);
            this.label2.TabIndex = 3;
            this.label2.Text = "target.csv File";
            // 
            // outputFolderBox
            // 
            this.outputFolderBox.Location = new System.Drawing.Point(6, 143);
            this.outputFolderBox.Name = "outputFolderBox";
            this.outputFolderBox.Size = new System.Drawing.Size(452, 20);
            this.outputFolderBox.TabIndex = 4;
            // 
            // label3
            // 
            this.label3.AutoSize = true;
            this.label3.Location = new System.Drawing.Point(6, 127);
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
            this.rtWindow.Location = new System.Drawing.Point(186, 336);
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
            this.label4.Location = new System.Drawing.Point(183, 320);
            this.label4.Name = "label4";
            this.label4.Size = new System.Drawing.Size(47, 13);
            this.label4.TabIndex = 7;
            this.label4.Text = "RT (min)";
            // 
            // RAWBrowse
            // 
            this.RAWBrowse.Location = new System.Drawing.Point(485, 90);
            this.RAWBrowse.Name = "RAWBrowse";
            this.RAWBrowse.Size = new System.Drawing.Size(75, 23);
            this.RAWBrowse.TabIndex = 12;
            this.RAWBrowse.Text = "Browse";
            this.RAWBrowse.UseVisualStyleBackColor = true;
            this.RAWBrowse.Click += new System.EventHandler(this.RAWBrowse_Click);
            // 
            // CSVBrowse
            // 
            this.CSVBrowse.Location = new System.Drawing.Point(485, 21);
            this.CSVBrowse.Name = "CSVBrowse";
            this.CSVBrowse.Size = new System.Drawing.Size(75, 23);
            this.CSVBrowse.TabIndex = 13;
            this.CSVBrowse.Text = "Browse";
            this.CSVBrowse.UseVisualStyleBackColor = true;
            this.CSVBrowse.Click += new System.EventHandler(this.CSVBrowse_Click);
            // 
            // OutputBrowse
            // 
            this.OutputBrowse.Location = new System.Drawing.Point(485, 141);
            this.OutputBrowse.Name = "OutputBrowse";
            this.OutputBrowse.Size = new System.Drawing.Size(75, 23);
            this.OutputBrowse.TabIndex = 14;
            this.OutputBrowse.Text = "Browse";
            this.OutputBrowse.UseVisualStyleBackColor = true;
            this.OutputBrowse.Click += new System.EventHandler(this.OutputBrowse_Click);
            // 
            // start
            // 
            this.start.Location = new System.Drawing.Point(485, 206);
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
            this.browseTargetInput.Filter = "Target CSV |*.csv";
            this.browseTargetInput.Multiselect = true;
            this.browseTargetInput.FileOk += new System.ComponentModel.CancelEventHandler(this.browseTargetInput_FileOk);
            // 
            // noiseBandCap
            // 
            this.noiseBandCap.AutoSize = true;
            this.noiseBandCap.Location = new System.Drawing.Point(387, 333);
            this.noiseBandCap.Name = "noiseBandCap";
            this.noiseBandCap.Size = new System.Drawing.Size(183, 17);
            this.noiseBandCap.TabIndex = 18;
            this.noiseBandCap.Text = "Noise Band Cap Missing Channel";
            this.noiseBandCap.UseVisualStyleBackColor = true;
            // 
            // coalescence
            // 
            this.coalescence.AutoSize = true;
            this.coalescence.Location = new System.Drawing.Point(387, 356);
            this.coalescence.Name = "coalescence";
            this.coalescence.Size = new System.Drawing.Size(147, 17);
            this.coalescence.TabIndex = 19;
            this.coalescence.Text = "Track Peak Coalescence";
            this.coalescence.UseVisualStyleBackColor = true;
            // 
            // rawFileBox
            // 
            this.rawFileBox.Location = new System.Drawing.Point(6, 93);
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
            this.signalToNoiseThreshold.Location = new System.Drawing.Point(249, 336);
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
            this.label6.Location = new System.Drawing.Point(246, 320);
            this.label6.Name = "label6";
            this.label6.Size = new System.Drawing.Size(50, 13);
            this.label6.TabIndex = 22;
            this.label6.Text = "Min. S/N";
            // 
            // label5
            // 
            this.label5.AutoSize = true;
            this.label5.Location = new System.Drawing.Point(304, 320);
            this.label5.Name = "label5";
            this.label5.Size = new System.Drawing.Size(47, 13);
            this.label5.TabIndex = 23;
            this.label5.Text = "Isotopes";
            // 
            // Isotopes
            // 
            this.Isotopes.Location = new System.Drawing.Point(307, 338);
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
            this.Conversion.Location = new System.Drawing.Point(387, 310);
            this.Conversion.Name = "Conversion";
            this.Conversion.Size = new System.Drawing.Size(157, 17);
            this.Conversion.TabIndex = 25;
            this.Conversion.Text = "Check for Label Conversion";
            this.Conversion.UseVisualStyleBackColor = true;
            // 
            // label7
            // 
            this.label7.AutoSize = true;
            this.label7.Location = new System.Drawing.Point(384, 294);
            this.label7.Name = "label7";
            this.label7.Size = new System.Drawing.Size(111, 13);
            this.label7.TabIndex = 26;
            this.label7.Text = "Extra Analysis Options";
            // 
            // label8
            // 
            this.label8.AutoSize = true;
            this.label8.Location = new System.Drawing.Point(6, 177);
            this.label8.Name = "label8";
            this.label8.Size = new System.Drawing.Size(96, 13);
            this.label8.TabIndex = 28;
            this.label8.Text = "Metabolic Labeling";
            // 
            // label9
            // 
            this.label9.AutoSize = true;
            this.label9.Location = new System.Drawing.Point(213, 177);
            this.label9.Name = "label9";
            this.label9.Size = new System.Drawing.Size(93, 13);
            this.label9.TabIndex = 29;
            this.label9.Text = "Chemical Labeling";
            // 
            // SILACLys8CN
            // 
            this.SILACLys8CN.AutoSize = true;
            this.SILACLys8CN.Location = new System.Drawing.Point(9, 206);
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
            this.SILACLys8D.Location = new System.Drawing.Point(9, 229);
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
            this.NeuCodeLys1.Location = new System.Drawing.Point(9, 295);
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
            this.NeuCodeLys8Duplex.Location = new System.Drawing.Point(9, 318);
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
            this.NeuCodeLys8Triplex.Location = new System.Drawing.Point(9, 341);
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
            this.NeuCodeLys8Fourplex.Location = new System.Drawing.Point(9, 364);
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
            this.NeuCodeLys8Sixplex.Location = new System.Drawing.Point(9, 387);
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
            this.label10.Location = new System.Drawing.Point(384, 177);
            this.label10.Name = "label10";
            this.label10.Size = new System.Drawing.Size(73, 13);
            this.label10.TabIndex = 37;
            this.label10.Text = "Cluster Labels";
            // 
            // NeuCodeLeu7Duplex
            // 
            this.NeuCodeLeu7Duplex.AutoSize = true;
            this.NeuCodeLeu7Duplex.Location = new System.Drawing.Point(9, 410);
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
            this.SILACLeu7CN.Location = new System.Drawing.Point(9, 252);
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
            this.SILACLeu7D.Location = new System.Drawing.Point(9, 275);
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
            this.CarbamylCN.Location = new System.Drawing.Point(216, 202);
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
            this.FourplexL.Location = new System.Drawing.Point(216, 225);
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
            this.FourplexM.Location = new System.Drawing.Point(216, 248);
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
            this.FourplexH.Location = new System.Drawing.Point(216, 271);
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
            this.Twelveplex.Location = new System.Drawing.Point(216, 294);
            this.Twelveplex.Name = "Twelveplex";
            this.Twelveplex.Size = new System.Drawing.Size(56, 17);
            this.Twelveplex.TabIndex = 45;
            this.Twelveplex.TabStop = true;
            this.Twelveplex.Text = "12plex";
            this.Twelveplex.UseVisualStyleBackColor = true;
            // 
            // IncompleteIncorporation
            // 
            this.IncompleteIncorporation.AutoSize = true;
            this.IncompleteIncorporation.Location = new System.Drawing.Point(387, 380);
            this.IncompleteIncorporation.Name = "IncompleteIncorporation";
            this.IncompleteIncorporation.Size = new System.Drawing.Size(169, 17);
            this.IncompleteIncorporation.TabIndex = 49;
            this.IncompleteIncorporation.Text = "Check for Partial Incorporation";
            this.IncompleteIncorporation.UseVisualStyleBackColor = true;
            // 
            // label11
            // 
            this.label11.AutoSize = true;
            this.label11.Location = new System.Drawing.Point(183, 364);
            this.label11.Name = "label11";
            this.label11.Size = new System.Drawing.Size(38, 13);
            this.label11.TabIndex = 50;
            this.label11.Text = "FWxM";
            // 
            // label12
            // 
            this.label12.AutoSize = true;
            this.label12.Location = new System.Drawing.Point(246, 363);
            this.label12.Name = "label12";
            this.label12.Size = new System.Drawing.Size(73, 13);
            this.label12.TabIndex = 51;
            this.label12.Text = "Resolution (K)";
            // 
            // PeakSeparation
            // 
            this.PeakSeparation.Location = new System.Drawing.Point(186, 380);
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
            this.QuantResolution.Location = new System.Drawing.Point(251, 380);
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
            // mTRAQ
            // 
            this.mTRAQ.AutoSize = true;
            this.mTRAQ.Location = new System.Drawing.Point(390, 208);
            this.mTRAQ.Name = "mTRAQ";
            this.mTRAQ.Size = new System.Drawing.Size(64, 17);
            this.mTRAQ.TabIndex = 54;
            this.mTRAQ.Text = "mTRAQ";
            this.mTRAQ.UseVisualStyleBackColor = true;
            // 
            // Arg
            // 
            this.Arg.AutoSize = true;
            this.Arg.Location = new System.Drawing.Point(390, 231);
            this.Arg.Name = "Arg";
            this.Arg.Size = new System.Drawing.Size(42, 17);
            this.Arg.TabIndex = 55;
            this.Arg.Text = "Arg";
            this.Arg.UseVisualStyleBackColor = true;
            // 
            // Leu
            // 
            this.Leu.AutoSize = true;
            this.Leu.Location = new System.Drawing.Point(390, 254);
            this.Leu.Name = "Leu";
            this.Leu.Size = new System.Drawing.Size(44, 17);
            this.Leu.TabIndex = 56;
            this.Leu.Text = "Leu";
            this.Leu.UseVisualStyleBackColor = true;
            // 
            // label13
            // 
            this.label13.AutoSize = true;
            this.label13.Location = new System.Drawing.Point(183, 410);
            this.label13.Name = "label13";
            this.label13.Size = new System.Drawing.Size(122, 13);
            this.label13.TabIndex = 57;
            this.label13.Text = "Quant Scan Alternatives";
            // 
            // MultipleInjections
            // 
            this.MultipleInjections.AutoSize = true;
            this.MultipleInjections.Location = new System.Drawing.Point(186, 426);
            this.MultipleInjections.Name = "MultipleInjections";
            this.MultipleInjections.Size = new System.Drawing.Size(77, 17);
            this.MultipleInjections.TabIndex = 59;
            this.MultipleInjections.Text = "Multi Inject";
            this.MultipleInjections.UseVisualStyleBackColor = true;
            // 
            // AGCBins
            // 
            this.AGCBins.AutoSize = true;
            this.AGCBins.Location = new System.Drawing.Point(186, 449);
            this.AGCBins.Name = "AGCBins";
            this.AGCBins.Size = new System.Drawing.Size(86, 17);
            this.AGCBins.TabIndex = 60;
            this.AGCBins.Text = "AGC Binning";
            this.AGCBins.UseVisualStyleBackColor = true;
            // 
            // Icat
            // 
            this.Icat.AutoSize = true;
            this.Icat.Location = new System.Drawing.Point(387, 410);
            this.Icat.Name = "Icat";
            this.Icat.Size = new System.Drawing.Size(49, 17);
            this.Icat.TabIndex = 61;
            this.Icat.TabStop = true;
            this.Icat.Text = "ICAT";
            this.Icat.UseVisualStyleBackColor = true;
            // 
            // listBox1
            // 
            this.listBox1.FormattingEnabled = true;
            this.listBox1.Location = new System.Drawing.Point(9, 21);
            this.listBox1.Name = "listBox1";
            this.listBox1.Size = new System.Drawing.Size(449, 56);
            this.listBox1.TabIndex = 62;
            // 
            // splitContainer1
            // 
            this.splitContainer1.Dock = System.Windows.Forms.DockStyle.Fill;
            this.splitContainer1.Location = new System.Drawing.Point(0, 0);
            this.splitContainer1.Name = "splitContainer1";
            this.splitContainer1.Orientation = System.Windows.Forms.Orientation.Horizontal;
            // 
            // splitContainer1.Panel1
            // 
            this.splitContainer1.Panel1.Controls.Add(this.button1);
            this.splitContainer1.Panel1.Controls.Add(this.IncompleteIncorporation);
            this.splitContainer1.Panel1.Controls.Add(this.listBox1);
            this.splitContainer1.Panel1.Controls.Add(this.label1);
            this.splitContainer1.Panel1.Controls.Add(this.Icat);
            this.splitContainer1.Panel1.Controls.Add(this.label2);
            this.splitContainer1.Panel1.Controls.Add(this.AGCBins);
            this.splitContainer1.Panel1.Controls.Add(this.outputFolderBox);
            this.splitContainer1.Panel1.Controls.Add(this.MultipleInjections);
            this.splitContainer1.Panel1.Controls.Add(this.label3);
            this.splitContainer1.Panel1.Controls.Add(this.label13);
            this.splitContainer1.Panel1.Controls.Add(this.rtWindow);
            this.splitContainer1.Panel1.Controls.Add(this.Leu);
            this.splitContainer1.Panel1.Controls.Add(this.label4);
            this.splitContainer1.Panel1.Controls.Add(this.Arg);
            this.splitContainer1.Panel1.Controls.Add(this.RAWBrowse);
            this.splitContainer1.Panel1.Controls.Add(this.mTRAQ);
            this.splitContainer1.Panel1.Controls.Add(this.CSVBrowse);
            this.splitContainer1.Panel1.Controls.Add(this.QuantResolution);
            this.splitContainer1.Panel1.Controls.Add(this.OutputBrowse);
            this.splitContainer1.Panel1.Controls.Add(this.PeakSeparation);
            this.splitContainer1.Panel1.Controls.Add(this.start);
            this.splitContainer1.Panel1.Controls.Add(this.label12);
            this.splitContainer1.Panel1.Controls.Add(this.noiseBandCap);
            this.splitContainer1.Panel1.Controls.Add(this.label11);
            this.splitContainer1.Panel1.Controls.Add(this.coalescence);
            this.splitContainer1.Panel1.Controls.Add(this.rawFileBox);
            this.splitContainer1.Panel1.Controls.Add(this.Twelveplex);
            this.splitContainer1.Panel1.Controls.Add(this.signalToNoiseThreshold);
            this.splitContainer1.Panel1.Controls.Add(this.FourplexH);
            this.splitContainer1.Panel1.Controls.Add(this.label6);
            this.splitContainer1.Panel1.Controls.Add(this.FourplexM);
            this.splitContainer1.Panel1.Controls.Add(this.label5);
            this.splitContainer1.Panel1.Controls.Add(this.FourplexL);
            this.splitContainer1.Panel1.Controls.Add(this.Isotopes);
            this.splitContainer1.Panel1.Controls.Add(this.CarbamylCN);
            this.splitContainer1.Panel1.Controls.Add(this.Conversion);
            this.splitContainer1.Panel1.Controls.Add(this.SILACLeu7D);
            this.splitContainer1.Panel1.Controls.Add(this.label7);
            this.splitContainer1.Panel1.Controls.Add(this.SILACLeu7CN);
            this.splitContainer1.Panel1.Controls.Add(this.label8);
            this.splitContainer1.Panel1.Controls.Add(this.NeuCodeLeu7Duplex);
            this.splitContainer1.Panel1.Controls.Add(this.label9);
            this.splitContainer1.Panel1.Controls.Add(this.label10);
            this.splitContainer1.Panel1.Controls.Add(this.SILACLys8CN);
            this.splitContainer1.Panel1.Controls.Add(this.NeuCodeLys8Sixplex);
            this.splitContainer1.Panel1.Controls.Add(this.SILACLys8D);
            this.splitContainer1.Panel1.Controls.Add(this.NeuCodeLys8Fourplex);
            this.splitContainer1.Panel1.Controls.Add(this.NeuCodeLys1);
            this.splitContainer1.Panel1.Controls.Add(this.NeuCodeLys8Triplex);
            this.splitContainer1.Panel1.Controls.Add(this.NeuCodeLys8Duplex);
            // 
            // splitContainer1.Panel2
            // 
            this.splitContainer1.Panel2.Controls.Add(this.richTextBox1);
            this.splitContainer1.Size = new System.Drawing.Size(622, 605);
            this.splitContainer1.SplitterDistance = 470;
            this.splitContainer1.TabIndex = 63;
            // 
            // richTextBox1
            // 
            this.richTextBox1.Dock = System.Windows.Forms.DockStyle.Fill;
            this.richTextBox1.Location = new System.Drawing.Point(0, 0);
            this.richTextBox1.Name = "richTextBox1";
            this.richTextBox1.Size = new System.Drawing.Size(622, 131);
            this.richTextBox1.TabIndex = 0;
            this.richTextBox1.Text = "";
            // 
            // button1
            // 
            this.button1.Location = new System.Drawing.Point(485, 54);
            this.button1.Name = "button1";
            this.button1.Size = new System.Drawing.Size(75, 23);
            this.button1.TabIndex = 63;
            this.button1.Text = "Clear";
            this.button1.UseVisualStyleBackColor = true;
            this.button1.Click += new System.EventHandler(this.button1_Click);
            // 
            // Form1
            // 
            this.AllowDrop = true;
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 13F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(622, 605);
            this.Controls.Add(this.splitContainer1);
            this.Name = "Form1";
            this.Text = "NeuQuant";
            this.Load += new System.EventHandler(this.Form1_Load);
            this.DragDrop += new System.Windows.Forms.DragEventHandler(this.Form1_DragDrop);
            this.DragEnter += new System.Windows.Forms.DragEventHandler(this.Form1_DragEnter);
            ((System.ComponentModel.ISupportInitialize)(this.rtWindow)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.signalToNoiseThreshold)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.Isotopes)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.PeakSeparation)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.QuantResolution)).EndInit();
            this.splitContainer1.Panel1.ResumeLayout(false);
            this.splitContainer1.Panel1.PerformLayout();
            this.splitContainer1.Panel2.ResumeLayout(false);
            ((System.ComponentModel.ISupportInitialize)(this.splitContainer1)).EndInit();
            this.splitContainer1.ResumeLayout(false);
            this.ResumeLayout(false);

        }

        #endregion

        private System.Windows.Forms.Label label1;
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
        private System.Windows.Forms.CheckBox IncompleteIncorporation;
        private System.Windows.Forms.Label label11;
        private System.Windows.Forms.Label label12;
        private System.Windows.Forms.NumericUpDown PeakSeparation;
        private System.Windows.Forms.NumericUpDown QuantResolution;
        private System.Windows.Forms.CheckBox mTRAQ;
        private System.Windows.Forms.CheckBox Arg;
        private System.Windows.Forms.CheckBox Leu;
        private System.Windows.Forms.Label label13;
        private System.Windows.Forms.CheckBox MultipleInjections;
        private System.Windows.Forms.CheckBox AGCBins;
        private System.Windows.Forms.RadioButton Icat;
        private System.Windows.Forms.ListBox listBox1;
        private System.Windows.Forms.SplitContainer splitContainer1;
        private System.Windows.Forms.RichTextBox richTextBox1;
        private System.Windows.Forms.Button button1;
    }
}

