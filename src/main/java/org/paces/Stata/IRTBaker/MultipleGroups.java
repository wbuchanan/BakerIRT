package org.paces.Stata.IRTBaker;

import javax.swing.*;

/**
 * Created by billy on 9/27/15.
 */
public class MultipleGroups {


	/**
	 * AppendixI.java
	 * Latent Population parameter estimation.
	 */
	static final int NL = 32; // Number of theta Levels
	static final int NI = 5; // Number of items I
	static final int NQ = 20; // Number of Quadratures
	static final int NC = 50; // Number of em Cycles
	static final int NG = 0;
	static final int MIT = 50;
	static final int MNC = 5;

	public MultipleGroups(String[] args) {
		this.latentPopulationParameter();
		this.difMultiGroup();
		this.augmentedModel();
	}

	public void latentPopulationParameter() {
		double[] x = new double[NQ];
		double[] ax = new double[NQ];
		double[] hu = new double[NL];
		double[] p = new double[NI];
		double[] sump = new double[NQ];
		double[][] lux = new double[NL][NQ];
		double[][] pux = new double[NL][NQ];
		// item parameter values
		double[] a = new double[NI];
		a[0] = 1.0;
		a[1] = 1.0;
		a[2] = 1.0;
		a[3] = 1.0;
		a[4] = 1.0;
		double[] b = new double[NI];
		b[0] = -1.256;
		b[1] = 0.475;
		b[2] = 1.236;
		b[3] = 0.168;
		b[4] = -0.623;
		double sumr = 1000.0;
		double d = 0.5;
		// d = interval for Sheppard's correction
		double[][] u = {
				{0.0, 0.0, 0.0, 0.0, 0.0},
				{0.0, 0.0, 0.0, 0.0, 1.0},
				{0.0, 0.0, 0.0, 1.0, 0.0},
				{0.0, 0.0, 0.0, 1.0, 1.0},
				{0.0, 0.0, 1.0, 0.0, 0.0},
				{0.0, 0.0, 1.0, 0.0, 1.0},
				{0.0, 0.0, 1.0, 1.0, 0.0},
				{0.0, 0.0, 1.0, 1.0, 1.0},
				{0.0, 1.0, 0.0, 0.0, 0.0},
				{0.0, 1.0, 0.0, 0.0, 1.0},
				{0.0, 1.0, 0.0, 1.0, 0.0},
				{0.0, 1.0, 0.0, 1.0, 1.0},
				{0.0, 1.0, 1.0, 0.0, 0.0},
				{0.0, 1.0, 1.0, 0.0, 1.0},
				{0.0, 1.0, 1.0, 1.0, 0.0},
				{0.0, 1.0, 1.0, 1.0, 1.0},
				{1.0, 0.0, 0.0, 0.0, 0.0},
				{1.0, 0.0, 0.0, 0.0, 1.0},
				{1.0, 0.0, 0.0, 1.0, 0.0},
				{1.0, 0.0, 0.0, 1.0, 1.0},
				{1.0, 0.0, 1.0, 0.0, 0.0},
				{1.0, 0.0, 1.0, 0.0, 1.0},
				{1.0, 0.0, 1.0, 1.0, 0.0},
				{1.0, 0.0, 1.0, 1.0, 1.0},
				{1.0, 1.0, 0.0, 0.0, 0.0},
				{1.0, 1.0, 0.0, 0.0, 1.0},
				{1.0, 1.0, 0.0, 1.0, 0.0},
				{1.0, 1.0, 0.0, 1.0, 1.0},
				{1.0, 1.0, 1.0, 0.0, 0.0},
				{1.0, 1.0, 1.0, 0.0, 1.0},
				{1.0, 1.0, 1.0, 1.0, 0.0},
				{1.0, 1.0, 1.0, 1.0, 1.0}
		};
		String[] id = {
				"00000", "00001", "00010", "00011", "00100",
				"00101", "00110", "00111", "01000", "01001",
				"01010", "01011", "01100", "01101", "01110",
				"01111", "10000", "10001", "10010", "10011",
				"10100", "10101", "10110", "10111", "11000",
				"11001", "11010", "11011", "11100", "11101",
				"11110", "11111"
		};
		double[] r = {
				3, 6, 2, 11, 1, 1, 3, 4, 1, 8,
				0, 16, 0, 3, 2, 15, 10, 29, 14, 81,
				3, 28, 15, 80, 16, 56, 21, 173, 11, 61,
				28, 298
		};
		double mu0 = 0.0;
		double sigma0 = 1.0;
		for (int c = 0; c <= NC - 1; c++) {
			double qnormcon = 0.0;
			for (int q = 0; q <= NQ - 1; q++) {
				x[q] = -4.75 + (q - 1.0) * d;
				qnormcon += Math.exp(-(1.0 / 2.0) *
						Math.pow((x[q] - mu0) / sigma0, 2));
			} // end for (int q ...
			for (int q = 0; q <= NQ - 1; q++) {
				ax[q] = Math.exp(-(1.0 / 2.0) *
						Math.pow((x[q] - mu0) / sigma0, 2)) / qnormcon;
			} // end for (int q ...
			for (int l = 0; l <= NL - 1; l++) {
				hu[l] = 0.0;
				for (int q = 0; q <= NQ - 1; q++) {
					lux[l][q] = 1.0;
				} // end for (int q ...
			} // end for (int l ...
			for (int l = 0; l <= NL - 1; l++) {
				for (int q = 0; q <= NQ - 1; q++) {
					for (int i = 0; i <= NI - 1; i++) {
						p[i] = 1.0 / (1.0 + Math.exp(-a[i] * (x[q] - b[i])));
						if (u[l][i] == 0.0) p[i] = 1.0 - p[i];
						lux[l][q] *= p[i];
					} // end for (int i ...
				} // end for (int q ...
				hu[l] = 0.0;
				for (int q = 0; q <= NQ - 1; q++) {
					hu[l] += lux[l][q] * ax[q];
				} // end for (int q ...
			} // end for (int l ...
			for (int l = 0; l <= NL - 1; l++) {
				for (int q = 0; q <= NQ - 1; q++) {
					pux[l][q] = lux[l][q] * ax[q] / hu[l];
				} // end for (int q ...
			} // end for (int l ...
			double mu = 0.0;
			double var = 0.0;
			for (int q = 0; q <= NQ - 1; q++) {
				sump[q] = 0.0;
				for (int l = 0; l <= NL - 1; l++) {
					sump[q] += pux[l][q] * r[l];
				} // end for (int l ...
				mu += x[q] * sump[q];
			} // end for (int q ...
			mu /= sumr;
			for (int q = 0; q <= NQ - 1; q++) {
				var += Math.pow(x[q] - mu, 2) * sump[q];
			} // end for (int q ...
			var /= sumr;
			double sigma = Math.sqrt(var);
			System.out.println("Cycle= " + (c + 1) +
					"  MU= " + mu +
					"  VAR= " + var +
					"  SIGMA= " + sigma);
			mu0 = mu;
			sigma0 = sigma;
			if (c == NC - 1) {
				// Sheppard's correction for final estimates
				double shevar = var - Math.pow(d, 2) / 12.0;
				double shesigma = Math.sqrt(shevar);
				System.out.println("MU= " + mu +
						"  VAR= " + var +
						"  SIGMA= " + sigma);
				System.out.println("Sheppard VAR= " + shevar +
						"  Sheppard SIGMA= " + shesigma);
				// information matrix
				double summu2 = 0.0;
				double summuvar = 0.0;
				double sumvar2 = 0.0;
				for (int l = 0; l <= NL - 1; l++) {
					double tempmu2 = 0.0;
					double tempvar2 = 0.0;
					for (int q = 0; q <= NQ - 1; q++) {
						tempmu2 += pux[l][q] * ((x[q] - mu) / var);
						tempvar2 += pux[l][q] *
								(Math.pow(x[q] - mu, 2) - var) / (2.0 * Math.pow(var, 2));
					} // end for (int q ...
					summu2 += r[l] * tempmu2 * tempmu2;
					summuvar += r[l] * tempmu2 * tempvar2;
					sumvar2 += r[l] * tempvar2 * tempvar2;
				} // end for (int l ..
				System.out.println("I(MU)= " + summu2 +
						"  I(MU,VAR)= " + summuvar +
						"  I(VAR)= " + sumvar2);
				double det = summu2 * sumvar2 - summuvar * summuvar;
				double varmu = sumvar2 / det;
				double varvar = summu2 / det;
				double covmuvar = -summuvar / det;
				double semu = Math.sqrt(varmu);
				double sevar = Math.sqrt(varvar);
				System.out.println("SE(MU)= " + semu +
						"  SE(VAR)= " + sevar +
						"  COVAR= " + covmuvar);
			} else {
				// normal processing
			} // end if (c ... else
		} // end for (int c ...
	}

	/**
	 * AppendixI1.java
	 * The compact DIF model.
	 */

	public void difMultiGroup() {
		double[] sax = new double[NG];
		double[] sn = new double[NG];
		double[] muh = new double[NG];
		double[] gof = new double[NG];
		double[][] ax = new double[NG][NQ];
		double[][] pb = new double[NG][NL];
		double[][] xb = new double[NG][NL];
		// y=indicator
		double[][][] lx = new double[NG][NL][NQ];
		double[][][] rb = new double[NG][NI][NQ];
		double[][][] nb = new double[NG][NI][NQ];
		// initial values
		double[] zeta = new double[NI];
		zeta[0] = 0.0;
		zeta[1] = 0.0;
		zeta[2] = 0.0;
		zeta[3] = 0.0;
		double[] lambda = new double[NI];
		lambda[0] = 1.0;
		lambda[1] = 1.0;
		lambda[2] = 1.0;
		lambda[3] = 1.0;
		double[][][] u = {
				{
						{0.0, 0.0, 0.0, 0.0},
						{0.0, 0.0, 0.0, 1.0},
						{1.0, 0.0, 0.0, 0.0},
						{1.0, 0.0, 0.0, 1.0},
						{0.0, 1.0, 0.0, 0.0},
						{0.0, 1.0, 0.0, 1.0},
						{0.0, 0.0, 1.0, 0.0},
						{0.0, 0.0, 1.0, 1.0},
						{1.0, 1.0, 0.0, 0.0},
						{1.0, 1.0, 0.0, 1.0},
						{1.0, 0.0, 1.0, 0.0},
						{1.0, 0.0, 1.0, 1.0},
						{0.0, 1.0, 1.0, 0.0},
						{0.0, 1.0, 1.0, 1.0},
						{1.0, 1.0, 1.0, 0.0},
						{1.0, 1.0, 1.0, 1.0},
				},
				{
						{0.0, 0.0, 0.0, 0.0},
						{0.0, 0.0, 0.0, 1.0},
						{1.0, 0.0, 0.0, 0.0},
						{1.0, 0.0, 0.0, 1.0},
						{0.0, 1.0, 0.0, 0.0},
						{0.0, 1.0, 0.0, 1.0},
						{0.0, 0.0, 1.0, 0.0},
						{0.0, 0.0, 1.0, 1.0},
						{1.0, 1.0, 0.0, 0.0},
						{1.0, 1.0, 0.0, 1.0},
						{1.0, 0.0, 1.0, 0.0},
						{1.0, 0.0, 1.0, 1.0},
						{0.0, 1.0, 1.0, 0.0},
						{0.0, 1.0, 1.0, 1.0},
						{1.0, 1.0, 1.0, 0.0},
						{1.0, 1.0, 1.0, 1.0},
				}
		};
		double[][][] y = {
				{
						{1.0, 1.0, 1.0, 1.0},
						{1.0, 1.0, 1.0, 1.0},
						{1.0, 1.0, 1.0, 1.0},
						{1.0, 1.0, 1.0, 1.0},
						{1.0, 1.0, 1.0, 1.0},
						{1.0, 1.0, 1.0, 1.0},
						{1.0, 1.0, 1.0, 1.0},
						{1.0, 1.0, 1.0, 1.0},
						{1.0, 1.0, 1.0, 1.0},
						{1.0, 1.0, 1.0, 1.0},
						{1.0, 1.0, 1.0, 1.0},
						{1.0, 1.0, 1.0, 1.0},
						{1.0, 1.0, 1.0, 1.0},
						{1.0, 1.0, 1.0, 1.0},
						{1.0, 1.0, 1.0, 1.0},
						{1.0, 1.0, 1.0, 1.0},
				},
				{
						{1.0, 1.0, 1.0, 1.0},
						{1.0, 1.0, 1.0, 1.0},
						{1.0, 1.0, 1.0, 1.0},
						{1.0, 1.0, 1.0, 1.0},
						{1.0, 1.0, 1.0, 1.0},
						{1.0, 1.0, 1.0, 1.0},
						{1.0, 1.0, 1.0, 1.0},
						{1.0, 1.0, 1.0, 1.0},
						{1.0, 1.0, 1.0, 1.0},
						{1.0, 1.0, 1.0, 1.0},
						{1.0, 1.0, 1.0, 1.0},
						{1.0, 1.0, 1.0, 1.0},
						{1.0, 1.0, 1.0, 1.0},
						{1.0, 1.0, 1.0, 1.0},
						{1.0, 1.0, 1.0, 1.0},
						{1.0, 1.0, 1.0, 1.0},
				}
		};
		double[][] x = {
				{-4.5, -3.5, -2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 3.5, 4.5},
				{-4.5, -3.5, -2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 3.5, 4.5}
		};
		double[] mu = new double[NG];
		mu[0] = 0.0;
		mu[1] = 0.0;
		double[] s2 = new double[NG];
		s2[0] = 1.0;
		s2[1] = 1.0;
		double[][] r = {
				{29.0, 7.0, 50.0, 30.0, 15.0, 4.0, 6.0, 0.0, 67.0, 63.0,
						12.0, 10.0, 2.0, 6.0, 22.0, 51.0},
				{22.0, 10.0, 30.0, 27.0, 13.0, 14.0, 1.0, 1.0, 24.0, 54.0,
						5.0, 8.0, 1.0, 8.0, 10.0, 57.0}
		};
		int trall = -1; // set traces
		String yn = JOptionPane.showInputDialog
				("Trace all? Y/N");
		if (yn.equalsIgnoreCase("y")) {
			trall = 0;
		} else {
			// normal processing
		} // end if (yn ... else
		for (int nc = 1; nc <= MNC; nc++) {
			System.out.println("Cycle= " + nc);
			// clear sums each iteration
			sax[0] = 0.0;
			sax[1] = 0.0;
			for (int g = 0; g <= NG - 1; g++) {
				for (int q = 0; q <= NQ - 1; q++) {
					ax[g][q] = (1.0 / Math.sqrt(2.0 * Math.PI)) *
							Math.exp(-(1.0 / 2.0) *
									Math.pow(x[g][q] - mu[g], 2) / s2[g]);
					sax[g] += ax[g][q];
				} // end for (int q ...
			} // end for (int g ...
			for (int g = 0; g <= NG - 1; g++) {
				for (int q = 0; q <= NQ - 1; q++) {
					ax[g][q] /= sax[g];
				} // end for (int q ...
			} // end for (int g ...
			if (trall == 0) {
				for (int g = 0; g <= NG - 1; g++) {
					for (int q = 0; q <= NQ - 1; q++) {
						System.out.println("G= " + (g + 1) +
								"  Q= " + (q + 1) +
								"  X= " + x[g][q] +
								"  AX= " + ax[g][q]);
					} // end for (int q ...
				} // end for (int g ...
			} else {
				// normal processing
			} // end if (trall ... else
			for (int g = 0; g <= NG - 1; g++) {
				for (int l = 0; l <= NL - 1; l++) {
					pb[g][l] = 0.0;
					for (int q = 0; q <= NQ - 1; q++) {
						lx[g][l][q] = 1.0;
					} // end for (int q ...
				} // end for (int l ...
			} // end for (int g ...
			// e step
			for (int g = 0; g <= NG - 1; g++) {
				for (int l = 0; l <= NL - 1; l++) {
					for (int i = 0; i <= NI - 1; i++) {
						for (int q = 0; q <= NQ - 1; q++) {
							if (y[g][l][i] != 0.0) {
								double pq = 1.0 /
										(1.0 + Math.exp(-(zeta[i] + lambda[i] * x[g][q])));
								if (l == 0 && trall == 0) {
									System.out.println("G= " + (g + 1) +
											"  I= " + (i + 1) +
											"  Q= " + (q + 1) +
											"  PQ= " + pq);
								} else {
									// normal processing
								} // end if (trall ... else
								if (u[g][l][i] == 0.0) pq = 1.0 - pq;
								lx[g][l][q] *= pq;
							} else {
								// normal process
							} // end if (y[g][l][i] ... else
						} // end for (int q ...
					} // end for (int i ...
					pb[g][l] = 0.0;
					for (int q = 0; q <= NQ - 1; q++) {
						pb[g][l] += lx[g][l][q] * ax[g][q];
					} // end for (int q ...
				} // end for (int l ...
			} // end for (int g ...
			if (trall == 0) {
				for (int g = 0; g <= NG - 1; g++) {
					for (int l = 0; l <= NL - 1; l++) {
						System.out.println("G= " + (g + 1) +
								"  L= " + (l + 1) +
								"  PB= " + pb[g][l]);
					} // end for (int l ...
				} // end for (int g ...
			} else {
				// normal processing
			} // end if (trall ... else
			if (trall == 0) {
				for (int g = 0; g <= NG - 1; g++) {
					for (int l = 0; l <= NL - 1; l++) {
						for (int q = 0; q <= NQ - 1; q++) {
							System.out.println("G= " + (g + 1) +
									"  L= " + (l + 1) +
									"  Q= " + (q + 1) +
									"  LX= " + lx[g][l][q]);
						} // end for (int q ...
					} // end for (int l ...
				} // end for (int g ...
			} else {
				// normal processing
			} // end if (trall ... else
			// n bar and r bar
			for (int g = 0; g <= NG - 1; g++) {
				for (int i = 0; i <= NI - 1; i++) {
					for (int q = 0; q <= NQ - 1; q++) {
						rb[g][i][q] = 0.0;
						nb[g][i][q] = 0.0;
						for (int l = 0; l <= NL - 1; l++) {
							if (y[g][l][i] != 0.0) {
								rb[g][i][q] += u[g][l][i] * r[g][l] * lx[g][l][q]
										* ax[g][q] / pb[g][l];
								nb[g][i][q] += r[g][l] * lx[g][l][q] * ax[g][q]
										/ pb[g][l];
							} else {
								// normal processing
							} // end if (y[g][l][i] ... else
						} // end for (int l ...
					} // end for (int q ...
				} // end for (int i ...
			} // end for (int g ...
			if (trall == 0) {
				for (int g = 0; g <= NG - 1; g++) {
					for (int i = 0; i <= NI - 1; i++) {
						for (int q = 0; q <= NQ - 1; q++) {
							System.out.println("G= " + (g + 1) +
									"  I= " + (i + 1) +
									"  Q= " + (q + 1) +
									"  RB= " + rb[g][i][q] +
									"  NB= " + nb[g][i][q]);
						} // end for (int q ...
					} // end for (int i ...
				} // end for (int g ...
			} else {
				// normal processing
			} // end if (trall ... else
			// probit stage m step
			for (int i = 0; i <= NI - 1; i++) {
				for (int it = 1; it <= MIT; it++) {
					double ssz = 0.0;
					double ssl = 0.0;
					double sszz = 0.0;
					double sszl = 0.0;
					double ssll = 0.0;
					for (int g = 0; g <= NG - 1; g++) {
						for (int q = 0; q <= NQ - 1; q++) {
							if (nb[g][i][q] != 0.0) {
								double px = 1.0 / (1.0 +
										Math.exp(-(zeta[i] + lambda[i] * x[g][q])));
								double w = px * (1.0 - px);
								ssz += rb[g][i][q] - nb[g][i][q] * px;
								ssl += (rb[g][i][q] - nb[g][i][q] * px) * x[g][q];
								sszz -= nb[g][i][q] * w;
								sszl -= nb[g][i][q] * w * x[g][q];
								ssll -= nb[g][i][q] * w * x[g][q] * x[g][q];
							} else {
								// normal processing
							} // end if (nb[g][i][q] ... else
						} // end for (int q ...
					} // end for (int g ...
					double det = sszz * ssll - sszl * sszl;
					double dzeta = (ssll * ssz - sszl * ssl) / det;
					double dlambda = (-sszl * ssz + sszz * ssl) / det;
					zeta[i] -= dzeta;
					lambda[i] -= dlambda;
					if (trall == 0) {
						System.out.println("I= " + (i + 1) +
								"  IT= " + it +
								"  ZETA= " + zeta[i] +
								"  LAMBDA= " + lambda[i]);
					} else {
						// normal processing
					} // end if (trall ... else
					if (Math.abs(dzeta) <= 0.01 && Math.abs(dlambda) <= 0.01) {
						break; // exit for (int it ...
					} else {
						// normal processing
					} // end if (Math.abs(dzeta) ... else
				} // end for (int it ...
			} // end for (int i ...
			for (int i = 0; i <= NI - 1; i++) {
				double b = -zeta[i] / lambda[i];
				System.out.println("I= " + (i + 1) +
						"  ZETA= " + zeta[i] +
						"  LAMBDA= " + lambda[i] +
						"  B= " + b);
			} // end for (int i ...
			for (int g = 0; g <= NG - 1; g++) {
				sn[g] = 0.0;
				for (int l = 0; l <= NL - 1; l++) {
					sn[g] += r[g][l];
				} // end for (int l ...
			} // end for (int g ...
			// estimate muh[0] only
			muh[0] = 0.0;
			for (int l = 0; l <= NL - 1; l++) {
				xb[0][l] = 0.0;
				for (int q = 0; q <= NQ - 1; q++) {
					xb[0][l] += x[0][q] * lx[0][l][q] * ax[0][q];
				} // end for (int q ...
				xb[0][l] /= pb[0][l];
				muh[0] += r[0][l] * xb[0][l];
			} // end for (int l ...
			muh[0] /= sn[0];
			System.out.println("MUH(1)= " + muh[0]);
			System.out.println("I(1)= " + sn[0]);
			mu[0] = muh[0];
			// above do not update
		} // end for (int nc ...
		// likelihood ratio goodness-of-fit statistic
		for (int g = 0; g <= NG - 1; g++) {
			for (int l = 0; l <= NL - 1; l++) {
				for (int q = 0; q <= NQ - 1; q++) {
					lx[g][l][q] = 1.0;
				} // end for (int q ...
			} // end for (int l ...
		} // end for (int g ...
		for (int g = 0; g <= NG - 1; g++) {
			for (int l = 0; l <= NL - 1; l++) {
				for (int i = 0; i <= NI - 1; i++) {
					for (int q = 0; q <= NQ - 1; q++) {
						if (y[g][l][i] != 0.0) {
							double pq = 1.0 /
									(1.0 + Math.exp(-(zeta[i] + lambda[i] * x[g][q])));
							if (u[g][l][i] == 0.0) pq = 1.0 - pq;
							lx[g][l][q] *= pq;
						} else {
							// normal processing
						} // end if (y[g][l][i] ... else
					} // end for (int q ...
				} // end for (int i ...
				pb[g][l] = 0.0;
				for (int q = 0; q <= NQ - 1; q++) {
					pb[g][l] += lx[g][l][q] * ax[g][q];
				} // end for (int q ...
			} // end for (int l ...
		} // end for (int g ...
		for (int g = 0; g <= NG - 1; g++) {
			for (int l = 0; l <= NL - 1; l++) {
				System.out.println("G= " + (g + 1) +
						"  L= " + (l + 1) +
						"  R= " + r[g][l] +
						"  N= " + sn[g] +
						"  PB= " + pb[g][l] +
						"  I*PB= " + sn[g] * pb[g][l]);
			} // end for (int l ...
		} // end for (int g ...
		gof[0] = 0.0;
		gof[1] = 0.0;
		for (int g = 0; g <= NG - 1; g++) {
			for (int l = 0; l <= NL - 1; l++) {
				if (r[g][l] != 0.0) {
					gof[g] += 2.0 * r[g][l] *
							Math.log(r[g][l] / (sn[g] * pb[g][l]));
				} else {
					// normal processing
				} // end if (r[g][l] ... else
			} // end for (int l ...
		} // end for (int g ...
		for (int g = 0; g <= NG - 1; g++) {
			System.out.println("G= " + (g + 1) +
					"  GOF= " + gof[g]);
		} // end for (int g ...
		double lrgof = 0.0;
		for (int g = 0; g <= NG - 1; g++) {
			lrgof += gof[g];
		} // end for (int g ...
		System.out.println("LR Goodness-of-Fit= " + lrgof);
		System.exit(0); // last line
	} // end public static ...

	/**
	 * AppendixI2.java
	 * An aumented model.
	 */

	public void augmentedModel() {
		double[] sax = new double[NG];
		double[] sn = new double[NG];
		double[] muh = new double[NG];
		double[] gof = new double[NG];
		double[][] ax = new double[NG][NQ];
		double[][] pb = new double[NG][NL];
		double[][] xb = new double[NG][NL];
		// y=indicator
		double[][][] lx = new double[NG][NL][NQ];
		double[][][] rb = new double[NG][NI][NQ];
		double[][][] nb = new double[NG][NI][NQ];
		// initial values
		double[] zeta = new double[NI];
		zeta[0] = 0.0;
		zeta[1] = 0.0;
		zeta[2] = 0.0;
		zeta[3] = 0.0;
		zeta[4] = 0.0;
		double[] lambda = new double[NI];
		lambda[0] = 1.0;
		lambda[1] = 1.0;
		lambda[2] = 1.0;
		lambda[3] = 1.0;
		lambda[4] = 1.0;
		double[][][] u = {
				{
						{0.0, 0.0, 0.0, 0.0, 0.0},
						{0.0, 0.0, 0.0, 0.0, 1.0},
						{1.0, 0.0, 0.0, 0.0, 0.0},
						{1.0, 0.0, 0.0, 0.0, 1.0},
						{0.0, 1.0, 0.0, 0.0, 0.0},
						{0.0, 1.0, 0.0, 0.0, 1.0},
						{0.0, 0.0, 1.0, 0.0, 0.0},
						{0.0, 0.0, 1.0, 0.0, 1.0},
						{1.0, 1.0, 0.0, 0.0, 0.0},
						{1.0, 1.0, 0.0, 0.0, 1.0},
						{1.0, 0.0, 1.0, 0.0, 0.0},
						{1.0, 0.0, 1.0, 0.0, 1.0},
						{0.0, 1.0, 1.0, 0.0, 0.0},
						{0.0, 1.0, 1.0, 0.0, 1.0},
						{1.0, 1.0, 1.0, 0.0, 0.0},
						{1.0, 1.0, 1.0, 0.0, 1.0},
				},
				{
						{0.0, 0.0, 0.0, 0.0, 0.0},
						{0.0, 0.0, 0.0, 1.0, 0.0},
						{1.0, 0.0, 0.0, 0.0, 0.0},
						{1.0, 0.0, 0.0, 1.0, 0.0},
						{0.0, 1.0, 0.0, 0.0, 0.0},
						{0.0, 1.0, 0.0, 1.0, 0.0},
						{0.0, 0.0, 1.0, 0.0, 0.0},
						{0.0, 0.0, 1.0, 1.0, 0.0},
						{1.0, 1.0, 0.0, 0.0, 0.0},
						{1.0, 1.0, 0.0, 1.0, 0.0},
						{1.0, 0.0, 1.0, 0.0, 0.0},
						{1.0, 0.0, 1.0, 1.0, 0.0},
						{0.0, 1.0, 1.0, 0.0, 0.0},
						{0.0, 1.0, 1.0, 1.0, 0.0},
						{1.0, 1.0, 1.0, 0.0, 0.0},
						{1.0, 1.0, 1.0, 1.0, 0.0},
				}
		};
		double[][][] y = {
				{
						{1.0, 1.0, 1.0, 0.0, 1.0},
						{1.0, 1.0, 1.0, 0.0, 1.0},
						{1.0, 1.0, 1.0, 0.0, 1.0},
						{1.0, 1.0, 1.0, 0.0, 1.0},
						{1.0, 1.0, 1.0, 0.0, 1.0},
						{1.0, 1.0, 1.0, 0.0, 1.0},
						{1.0, 1.0, 1.0, 0.0, 1.0},
						{1.0, 1.0, 1.0, 0.0, 1.0},
						{1.0, 1.0, 1.0, 0.0, 1.0},
						{1.0, 1.0, 1.0, 0.0, 1.0},
						{1.0, 1.0, 1.0, 0.0, 1.0},
						{1.0, 1.0, 1.0, 0.0, 1.0},
						{1.0, 1.0, 1.0, 0.0, 1.0},
						{1.0, 1.0, 1.0, 0.0, 1.0},
						{1.0, 1.0, 1.0, 0.0, 1.0},
						{1.0, 1.0, 1.0, 0.0, 1.0},
				},
				{
						{1.0, 1.0, 1.0, 1.0, 0.0},
						{1.0, 1.0, 1.0, 1.0, 0.0},
						{1.0, 1.0, 1.0, 1.0, 0.0},
						{1.0, 1.0, 1.0, 1.0, 0.0},
						{1.0, 1.0, 1.0, 1.0, 0.0},
						{1.0, 1.0, 1.0, 1.0, 0.0},
						{1.0, 1.0, 1.0, 1.0, 0.0},
						{1.0, 1.0, 1.0, 1.0, 0.0},
						{1.0, 1.0, 1.0, 1.0, 0.0},
						{1.0, 1.0, 1.0, 1.0, 0.0},
						{1.0, 1.0, 1.0, 1.0, 0.0},
						{1.0, 1.0, 1.0, 1.0, 0.0},
						{1.0, 1.0, 1.0, 1.0, 0.0},
						{1.0, 1.0, 1.0, 1.0, 0.0},
						{1.0, 1.0, 1.0, 1.0, 0.0},
						{1.0, 1.0, 1.0, 1.0, 0.0},
				}
		};
		double[][] x = {
				{-4.5, -3.5, -2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 3.5, 4.5},
				{-4.5, -3.5, -2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 3.5, 4.5}
		};
		double[] mu = new double[NG];
		mu[0] = 0.0;
		mu[1] = 0.0;
		double[] s2 = new double[NG];
		s2[0] = 1.0;
		s2[1] = 1.0;
		double[][] r = {
				{29.0, 7.0, 50.0, 30.0, 15.0, 4.0, 6.0, 0.0, 67.0, 63.0,
						12.0, 10.0, 2.0, 6.0, 22.0, 51.0},
				{22.0, 10.0, 30.0, 27.0, 13.0, 14.0, 1.0, 1.0, 24.0, 54.0,
						5.0, 8.0, 1.0, 8.0, 10.0, 57.0}
		};
		int trall = -1; // set traces
		String yn = JOptionPane.showInputDialog
				("Trace all? Y/N");
		if (yn.equalsIgnoreCase("y")) {
			trall = 0;
		} else {
			// normal processing
		} // end if (yn ... else
		for (int nc = 1; nc <= MNC; nc++) {
			System.out.println("Cycle= " + nc);
			// clear sums each iteration
			sax[0] = 0.0;
			sax[1] = 0.0;
			for (int g = 0; g <= NG - 1; g++) {
				for (int q = 0; q <= NQ - 1; q++) {
					ax[g][q] = (1.0 / Math.sqrt(2.0 * Math.PI)) *
							Math.exp(-(1.0 / 2.0) *
									Math.pow(x[g][q] - mu[g], 2) / s2[g]);
					sax[g] += ax[g][q];
				} // end for (int q ...
			} // end for (int g ...
			for (int g = 0; g <= NG - 1; g++) {
				for (int q = 0; q <= NQ - 1; q++) {
					ax[g][q] /= sax[g];
				} // end for (int q ...
			} // end for (int g ...
			if (trall == 0) {
				for (int g = 0; g <= NG - 1; g++) {
					for (int q = 0; q <= NQ - 1; q++) {
						System.out.println("G= " + (g + 1) +
								"  Q= " + (q + 1) +
								"  X= " + x[g][q] +
								"  AX= " + ax[g][q]);
					} // end for (int q ...
				} // end for (int g ...
			} else {
				// normal processing
			} // end if (trall ... else
			for (int g = 0; g <= NG - 1; g++) {
				for (int l = 0; l <= NL - 1; l++) {
					pb[g][l] = 0.0;
					for (int q = 0; q <= NQ - 1; q++) {
						lx[g][l][q] = 1.0;
					} // end for (int q ...
				} // end for (int l ...
			} // end for (int g ...
			// e step
			for (int g = 0; g <= NG - 1; g++) {
				for (int l = 0; l <= NL - 1; l++) {
					for (int i = 0; i <= NI - 1; i++) {
						for (int q = 0; q <= NQ - 1; q++) {
							if (y[g][l][i] != 0.0) {
								double pq = 1.0 /
										(1.0 + Math.exp(-(zeta[i] + lambda[i] * x[g][q])));
								if (l == 0 && trall == 0) {
									System.out.println("G= " + (g + 1) +
											"  I= " + (i + 1) +
											"  Q= " + (q + 1) +
											"  PQ= " + pq);
								} else {
									// normal processing
								} // end if (trall ... else
								if (u[g][l][i] == 0.0) pq = 1.0 - pq;
								lx[g][l][q] *= pq;
							} else {
								// normal process
							} // end if (y[g][l][i] ... else
						} // end for (int q ...
					} // end for (int i ...
					pb[g][l] = 0.0;
					for (int q = 0; q <= NQ - 1; q++) {
						pb[g][l] += lx[g][l][q] * ax[g][q];
					} // end for (int q ...
				} // end for (int l ...
			} // end for (int g ...
			if (trall == 0) {
				for (int g = 0; g <= NG - 1; g++) {
					for (int l = 0; l <= NL - 1; l++) {
						System.out.println("G= " + (g + 1) +
								"  L= " + (l + 1) +
								"  PB= " + pb[g][l]);
					} // end for (int l ...
				} // end for (int g ...
			} else {
				// normal processing
			} // end if (trall ... else
			if (trall == 0) {
				for (int g = 0; g <= NG - 1; g++) {
					for (int l = 0; l <= NL - 1; l++) {
						for (int q = 0; q <= NQ - 1; q++) {
							System.out.println("G= " + (g + 1) +
									"  L= " + (l + 1) +
									"  Q= " + (q + 1) +
									"  LX= " + lx[g][l][q]);
						} // end for (int q ...
					} // end for (int l ...
				} // end for (int g ...
			} else {
				// normal processing
			} // end if (trall ... else
			// n bar and r bar
			for (int g = 0; g <= NG - 1; g++) {
				for (int i = 0; i <= NI - 1; i++) {
					for (int q = 0; q <= NQ - 1; q++) {
						rb[g][i][q] = 0.0;
						nb[g][i][q] = 0.0;
						for (int l = 0; l <= NL - 1; l++) {
							if (y[g][l][i] != 0.0) {
								rb[g][i][q] += u[g][l][i] * r[g][l] * lx[g][l][q]
										* ax[g][q] / pb[g][l];
								nb[g][i][q] += r[g][l] * lx[g][l][q] * ax[g][q]
										/ pb[g][l];
							} else {
								// normal processing
							} // end if (y[g][l][i] ... else
						} // end for (int l ...
					} // end for (int q ...
				} // end for (int i ...
			} // end for (int g ...
			if (trall == 0) {
				for (int g = 0; g <= NG - 1; g++) {
					for (int i = 0; i <= NI - 1; i++) {
						for (int q = 0; q <= NQ - 1; q++) {
							System.out.println("G= " + (g + 1) +
									"  I= " + (i + 1) +
									"  Q= " + (q + 1) +
									"  RB= " + rb[g][i][q] +
									"  NB= " + nb[g][i][q]);
						} // end for (int q ...
					} // end for (int i ...
				} // end for (int g ...
			} else {
				// normal processing
			} // end if (trall ... else
			// probit stage m step
			for (int i = 0; i <= NI - 1; i++) {
				for (int it = 1; it <= MIT; it++) {
					double ssz = 0.0;
					double ssl = 0.0;
					double sszz = 0.0;
					double sszl = 0.0;
					double ssll = 0.0;
					for (int g = 0; g <= NG - 1; g++) {
						for (int q = 0; q <= NQ - 1; q++) {
							if (nb[g][i][q] != 0.0) {
								double px = 1.0 / (1.0 +
										Math.exp(-(zeta[i] + lambda[i] * x[g][q])));
								double w = px * (1.0 - px);
								ssz += rb[g][i][q] - nb[g][i][q] * px;
								ssl += (rb[g][i][q] - nb[g][i][q] * px) * x[g][q];
								sszz -= nb[g][i][q] * w;
								sszl -= nb[g][i][q] * w * x[g][q];
								ssll -= nb[g][i][q] * w * x[g][q] * x[g][q];
							} else {
								// normal processing
							} // end if (nb[g][i][q] ... else
						} // end for (int q ...
					} // end for (int g ...
					double det = sszz * ssll - sszl * sszl;
					double dzeta = (ssll * ssz - sszl * ssl) / det;
					double dlambda = (-sszl * ssz + sszz * ssl) / det;
					zeta[i] -= dzeta;
					lambda[i] -= dlambda;
					if (trall == 0) {
						System.out.println("I= " + (i + 1) +
								"  IT= " + it +
								"  ZETA= " + zeta[i] +
								"  LAMBDA= " + lambda[i]);
					} else {
						// normal processing
					} // end if (trall ... else
					if (Math.abs(dzeta) <= 0.01 && Math.abs(dlambda) <= 0.01) {
						break; // exit for (int it ...
					} else {
						// normal processing
					} // end if (Math.abs(dzeta) ... else
				} // end for (int it ...
			} // end for (int i ...
			for (int i = 0; i <= NI - 1; i++) {
				double b = -zeta[i] / lambda[i];
				System.out.println("I= " + (i + 1) +
						"  ZETA= " + zeta[i] +
						"  LAMBDA= " + lambda[i] +
						"  B= " + b);
			} // end for (int i ...
			for (int g = 0; g <= NG - 1; g++) {
				sn[g] = 0.0;
				for (int l = 0; l <= NL - 1; l++) {
					sn[g] += r[g][l];
				} // end for (int l ...
			} // end for (int g ...
			// estimate muh[0] only
			muh[0] = 0.0;
			for (int l = 0; l <= NL - 1; l++) {
				xb[0][l] = 0.0;
				for (int q = 0; q <= NQ - 1; q++) {
					xb[0][l] += x[0][q] * lx[0][l][q] * ax[0][q];
				} // end for (int q ...
				xb[0][l] /= pb[0][l];
				muh[0] += r[0][l] * xb[0][l];
			} // end for (int l ...
			muh[0] /= sn[0];
			System.out.println("MUH(1)= " + muh[0]);
			System.out.println("I(1)= " + sn[0]);
			mu[0] = muh[0];
			// above do not update
		} // end for (int nc ...
		// likelihood ratio goodness-of-fit statistic
		for (int g = 0; g <= NG - 1; g++) {
			for (int l = 0; l <= NL - 1; l++) {
				for (int q = 0; q <= NQ - 1; q++) {
					lx[g][l][q] = 1.0;
				} // end for (int q ...
			} // end for (int l ...
		} // end for (int g ...
		for (int g = 0; g <= NG - 1; g++) {
			for (int l = 0; l <= NL - 1; l++) {
				for (int i = 0; i <= NI - 1; i++) {
					for (int q = 0; q <= NQ - 1; q++) {
						if (y[g][l][i] != 0.0) {
							double pq = 1.0 /
									(1.0 + Math.exp(-(zeta[i] + lambda[i] * x[g][q])));
							if (u[g][l][i] == 0.0) pq = 1.0 - pq;
							lx[g][l][q] *= pq;
						} else {
							// normal processing
						} // end if (y[g][l][i] ... else
					} // end for (int q ...
				} // end for (int i ...
				pb[g][l] = 0.0;
				for (int q = 0; q <= NQ - 1; q++) {
					pb[g][l] += lx[g][l][q] * ax[g][q];
				} // end for (int q ...
			} // end for (int l ...
		} // end for (int g ...
		for (int g = 0; g <= NG - 1; g++) {
			for (int l = 0; l <= NL - 1; l++) {
				System.out.println("G= " + (g + 1) +
						"  L= " + (l + 1) +
						"  R= " + r[g][l] +
						"  N= " + sn[g] +
						"  PB= " + pb[g][l] +
						"  I*PB= " + sn[g] * pb[g][l]);
			} // end for (int l ...
		} // end for (int g ...
		gof[0] = 0.0;
		gof[1] = 0.0;
		for (int g = 0; g <= NG - 1; g++) {
			for (int l = 0; l <= NL - 1; l++) {
				if (r[g][l] != 0.0) {
					gof[g] += 2.0 * r[g][l] *
							Math.log(r[g][l] / (sn[g] * pb[g][l]));
				} else {
					// normal processing
				} // end if (r[g][l] ... else
			} // end for (int l ...
		} // end for (int g ...
		for (int g = 0; g <= NG - 1; g++) {
			System.out.println("G= " + (g + 1) +
					"  GOF= " + gof[g]);
		} // end for (int g ...
		double lrgof = 0.0;
		for (int g = 0; g <= NG - 1; g++) {
			lrgof += gof[g];
		} // end for (int g ...
		System.out.println("LR Goodness-of-Fit= " + lrgof);
		System.exit(0); // last line
	} // end public static ...

}