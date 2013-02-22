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
            this.SILAC = new System.Windows.Forms.RadioButton();
            this.HILAC = new System.Windows.Forms.RadioButton();
            this.channels = new System.Windows.Forms.NumericUpDown();
            this.label5 = new System.Windows.Forms.Label();
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
            this.dNLMaximum = new System.Windows.Forms.TextBox();
            this.label7 = new System.Windows.Forms.Label();
            ((System.ComponentModel.ISupportInitialize)(this.rtWindow)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.channels)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.signalToNoiseThreshold)).BeginInit();
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
            this.csvInputBox.Size = new System.Drawing.Size(297, 20);
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
            this.outputFolderBox.Size = new System.Drawing.Size(300, 20);
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
            this.rtWindow.Location = new System.Drawing.Point(15, 210);
            this.rtWindow.Maximum = new decimal(new int[] {
            5,
            0,
            0,
            0});
            this.rtWindow.Name = "rtWindow";
            this.rtWindow.Size = new System.Drawing.Size(120, 20);
            this.rtWindow.TabIndex = 6;
            this.rtWindow.Value = new decimal(new int[] {
            1,
            0,
            0,
            0});
            // 
            // label4
            // 
            this.label4.AutoSize = true;
            this.label4.Location = new System.Drawing.Point(12, 194);
            this.label4.Name = "label4";
            this.label4.Size = new System.Drawing.Size(89, 13);
            this.label4.TabIndex = 7;
            this.label4.Text = "RT Window (min)";
            // 
            // SILAC
            // 
            this.SILAC.AutoSize = true;
            this.SILAC.Location = new System.Drawing.Point(210, 253);
            this.SILAC.Name = "SILAC";
            this.SILAC.Size = new System.Drawing.Size(107, 17);
            this.SILAC.TabIndex = 8;
            this.SILAC.TabStop = true;
            this.SILAC.Text = "Traditional SILAC";
            this.SILAC.UseVisualStyleBackColor = true;
            // 
            // HILAC
            // 
            this.HILAC.AutoSize = true;
            this.HILAC.Location = new System.Drawing.Point(210, 280);
            this.HILAC.Name = "HILAC";
            this.HILAC.Size = new System.Drawing.Size(90, 17);
            this.HILAC.TabIndex = 9;
            this.HILAC.TabStop = true;
            this.HILAC.Text = "OMNE SILAC";
            this.HILAC.UseVisualStyleBackColor = true;
            // 
            // channels
            // 
            this.channels.Increment = new decimal(new int[] {
            2,
            0,
            0,
            0});
            this.channels.Location = new System.Drawing.Point(160, 210);
            this.channels.Maximum = new decimal(new int[] {
            12,
            0,
            0,
            0});
            this.channels.Minimum = new decimal(new int[] {
            2,
            0,
            0,
            0});
            this.channels.Name = "channels";
            this.channels.Size = new System.Drawing.Size(38, 20);
            this.channels.TabIndex = 10;
            this.channels.Value = new decimal(new int[] {
            2,
            0,
            0,
            0});
            // 
            // label5
            // 
            this.label5.AutoSize = true;
            this.label5.Location = new System.Drawing.Point(157, 194);
            this.label5.Name = "label5";
            this.label5.Size = new System.Drawing.Size(61, 13);
            this.label5.TabIndex = 11;
            this.label5.Text = "# Channels";
            // 
            // RAWBrowse
            // 
            this.RAWBrowse.Location = new System.Drawing.Point(329, 46);
            this.RAWBrowse.Name = "RAWBrowse";
            this.RAWBrowse.Size = new System.Drawing.Size(75, 23);
            this.RAWBrowse.TabIndex = 12;
            this.RAWBrowse.Text = "Browse";
            this.RAWBrowse.UseVisualStyleBackColor = true;
            this.RAWBrowse.Click += new System.EventHandler(this.RAWBrowse_Click);
            // 
            // CSVBrowse
            // 
            this.CSVBrowse.Location = new System.Drawing.Point(329, 90);
            this.CSVBrowse.Name = "CSVBrowse";
            this.CSVBrowse.Size = new System.Drawing.Size(75, 23);
            this.CSVBrowse.TabIndex = 13;
            this.CSVBrowse.Text = "Browse";
            this.CSVBrowse.UseVisualStyleBackColor = true;
            this.CSVBrowse.Click += new System.EventHandler(this.CSVBrowse_Click);
            // 
            // OutputBrowse
            // 
            this.OutputBrowse.Location = new System.Drawing.Point(329, 145);
            this.OutputBrowse.Name = "OutputBrowse";
            this.OutputBrowse.Size = new System.Drawing.Size(75, 23);
            this.OutputBrowse.TabIndex = 14;
            this.OutputBrowse.Text = "Browse";
            this.OutputBrowse.UseVisualStyleBackColor = true;
            this.OutputBrowse.Click += new System.EventHandler(this.OutputBrowse_Click);
            // 
            // start
            // 
            this.start.Location = new System.Drawing.Point(322, 277);
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
            this.noiseBandCap.Location = new System.Drawing.Point(15, 254);
            this.noiseBandCap.Name = "noiseBandCap";
            this.noiseBandCap.Size = new System.Drawing.Size(183, 17);
            this.noiseBandCap.TabIndex = 18;
            this.noiseBandCap.Text = "Noise Band Cap Missing Channel";
            this.noiseBandCap.UseVisualStyleBackColor = true;
            // 
            // coalescence
            // 
            this.coalescence.AutoSize = true;
            this.coalescence.Location = new System.Drawing.Point(15, 277);
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
            this.rawFileBox.Size = new System.Drawing.Size(302, 20);
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
            this.signalToNoiseThreshold.Location = new System.Drawing.Point(238, 210);
            this.signalToNoiseThreshold.Maximum = new decimal(new int[] {
            10,
            0,
            0,
            0});
            this.signalToNoiseThreshold.Name = "signalToNoiseThreshold";
            this.signalToNoiseThreshold.Size = new System.Drawing.Size(51, 20);
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
            this.label6.Location = new System.Drawing.Point(224, 194);
            this.label6.Name = "label6";
            this.label6.Size = new System.Drawing.Size(71, 13);
            this.label6.TabIndex = 22;
            this.label6.Text = "Minimum S/N";
            // 
            // dNLMaximum
            // 
            this.dNLMaximum.Location = new System.Drawing.Point(304, 210);
            this.dNLMaximum.Name = "dNLMaximum";
            this.dNLMaximum.Size = new System.Drawing.Size(100, 20);
            this.dNLMaximum.TabIndex = 23;
            // 
            // label7
            // 
            this.label7.AutoSize = true;
            this.label7.Location = new System.Drawing.Point(301, 194);
            this.label7.Name = "label7";
            this.label7.Size = new System.Drawing.Size(74, 13);
            this.label7.TabIndex = 24;
            this.label7.Text = "Maximum dNL";
            // 
            // Form1
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 13F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(409, 315);
            this.Controls.Add(this.label7);
            this.Controls.Add(this.dNLMaximum);
            this.Controls.Add(this.label6);
            this.Controls.Add(this.signalToNoiseThreshold);
            this.Controls.Add(this.rawFileBox);
            this.Controls.Add(this.coalescence);
            this.Controls.Add(this.noiseBandCap);
            this.Controls.Add(this.start);
            this.Controls.Add(this.OutputBrowse);
            this.Controls.Add(this.CSVBrowse);
            this.Controls.Add(this.RAWBrowse);
            this.Controls.Add(this.label5);
            this.Controls.Add(this.channels);
            this.Controls.Add(this.HILAC);
            this.Controls.Add(this.SILAC);
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
            ((System.ComponentModel.ISupportInitialize)(this.channels)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.signalToNoiseThreshold)).EndInit();
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
        private System.Windows.Forms.RadioButton SILAC;
        private System.Windows.Forms.RadioButton HILAC;
        private System.Windows.Forms.NumericUpDown channels;
        private System.Windows.Forms.Label label5;
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
        private System.Windows.Forms.TextBox dNLMaximum;
        private System.Windows.Forms.Label label7;
    }
}

