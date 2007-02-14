namespace distribution_explorer
{
    partial class distexSplash
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
          System.ComponentModel.ComponentResourceManager resources = new System.ComponentModel.ComponentResourceManager(typeof(distexSplash));
          this.labelApplicationTitle = new System.Windows.Forms.Label();
          this.labelApplicationVersion = new System.Windows.Forms.Label();
          this.labelApplicationCopyright = new System.Windows.Forms.Label();
          this.labelApplicationDescription = new System.Windows.Forms.Label();
          this.SuspendLayout();
          // 
          // labelApplicationTitle
          // 
          this.labelApplicationTitle.AutoSize = true;
          this.labelApplicationTitle.BackColor = System.Drawing.Color.White;
          this.labelApplicationTitle.Font = new System.Drawing.Font("Microsoft Sans Serif", 24F);
          this.labelApplicationTitle.ForeColor = System.Drawing.Color.Black;
          this.labelApplicationTitle.Location = new System.Drawing.Point(30, 126);
          this.labelApplicationTitle.Name = "labelApplicationTitle";
          this.labelApplicationTitle.Size = new System.Drawing.Size(376, 46);
          this.labelApplicationTitle.TabIndex = 0;
          this.labelApplicationTitle.Text = "labelApplicationTitle";
          // 
          // labelApplicationVersion
          // 
          this.labelApplicationVersion.AutoSize = true;
          this.labelApplicationVersion.BackColor = System.Drawing.Color.White;
          this.labelApplicationVersion.Font = new System.Drawing.Font("Microsoft Sans Serif", 12F, System.Drawing.FontStyle.Regular, System.Drawing.GraphicsUnit.Point, ((byte)(0)));
          this.labelApplicationVersion.ForeColor = System.Drawing.Color.Black;
          this.labelApplicationVersion.Location = new System.Drawing.Point(33, 267);
          this.labelApplicationVersion.Name = "labelApplicationVersion";
          this.labelApplicationVersion.Size = new System.Drawing.Size(216, 25);
          this.labelApplicationVersion.TabIndex = 1;
          this.labelApplicationVersion.Text = "labelApplicationVersion";
          this.labelApplicationVersion.Click += new System.EventHandler(this.labelApplicationVersion_Click);
          // 
          // labelApplicationCopyright
          // 
          this.labelApplicationCopyright.AutoSize = true;
          this.labelApplicationCopyright.BackColor = System.Drawing.Color.White;
          this.labelApplicationCopyright.Font = new System.Drawing.Font("Microsoft Sans Serif", 12F);
          this.labelApplicationCopyright.ForeColor = System.Drawing.Color.Navy;
          this.labelApplicationCopyright.Location = new System.Drawing.Point(33, 362);
          this.labelApplicationCopyright.Name = "labelApplicationCopyright";
          this.labelApplicationCopyright.Size = new System.Drawing.Size(233, 25);
          this.labelApplicationCopyright.TabIndex = 2;
          this.labelApplicationCopyright.Text = "labelApplicationCopyright";
          // 
          // labelApplicationDescription
          // 
          this.labelApplicationDescription.AutoSize = true;
          this.labelApplicationDescription.Font = new System.Drawing.Font("Microsoft Sans Serif", 12F);
          this.labelApplicationDescription.Location = new System.Drawing.Point(33, 465);
          this.labelApplicationDescription.Name = "labelApplicationDescription";
          this.labelApplicationDescription.Size = new System.Drawing.Size(246, 25);
          this.labelApplicationDescription.TabIndex = 3;
          this.labelApplicationDescription.Text = "labelApplicationDescription";
          // 
          // distexSplash
          // 
          this.AutoScaleDimensions = new System.Drawing.SizeF(8F, 18F);
          this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
          this.BackColor = System.Drawing.Color.White;
          this.BackgroundImage = global::distribution_explorer.Properties.Resources.ToolkitLogo;
          this.ClientSize = new System.Drawing.Size(800, 675);
          this.Controls.Add(this.labelApplicationDescription);
          this.Controls.Add(this.labelApplicationCopyright);
          this.Controls.Add(this.labelApplicationVersion);
          this.Controls.Add(this.labelApplicationTitle);
          this.Font = new System.Drawing.Font("Tahoma", 9F, System.Drawing.FontStyle.Regular, System.Drawing.GraphicsUnit.Point, ((byte)(0)));
          this.FormBorderStyle = System.Windows.Forms.FormBorderStyle.None;
          this.Icon = ((System.Drawing.Icon)(resources.GetObject("$this.Icon")));
          this.Name = "distexSplash";
          this.ShowInTaskbar = false;
          this.StartPosition = System.Windows.Forms.FormStartPosition.CenterScreen;
          this.Text = "Statistical Distribution Explorer";
          this.ResumeLayout(false);
          this.PerformLayout();

        }

        #endregion

        private System.Windows.Forms.Label labelApplicationTitle;
        private System.Windows.Forms.Label labelApplicationVersion;
      private System.Windows.Forms.Label labelApplicationCopyright;
      private System.Windows.Forms.Label labelApplicationDescription;
    }
}