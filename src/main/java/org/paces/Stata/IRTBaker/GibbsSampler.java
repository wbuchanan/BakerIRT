package org.paces.Stata.IRTBaker;

import javax.swing.*;
import java.util.Random;

/**
 * Created by billy on 9/27/15.
 */
public class GibbsSampler {

/**
 * AppendixK.java
 * Gibbs sampler for Memory data.
 */
		static final int NEXAM = 100; // Number of EXAMinees
		static final int NITEM = 50; // Number of ITEMs
		static final int NP = 2; // Number of Parameters
		public GibbsSampler(String[] args) {
			double[] a = new double[NITEM];
			double[] g = new double[NITEM];
			double[] th = new double[NEXAM];
			double[][] z = new double[NEXAM][NITEM];
			double[][] xpx = new double[NP][NP];
			double[][] xpxinv = new double[NP][NP];
			double[][] amat = new double[NP][NP];
			double[][] amatt = new double[NP][NP];
			double[][] pp = new double[NP][NP];
			String title = JOptionPane.showInputDialog
					("Enter title for run. ");
			System.out.println(title);
			String nexam = JOptionPane.showInputDialog
					("Enter number of examinees. ");
			int NI = Integer.parseInt(nexam);
			System.out.println("Number of examinees= " + NI);
			String nitem = JOptionPane.showInputDialog
					("Enter number of items. ");
			int NJ = Integer.parseInt(nitem);
			System.out.println("Number of items= " + NJ);
			double[][] y = {
					{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
					{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
					{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
					{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
					{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
					{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 },
					{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0 },
					{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0 },
					{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0 },
					{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0 },
					{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0 },
					{ 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 },
					{ 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0 },
					{ 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0 },
					{ 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 },
					{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0 },
					{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0 },
					{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0 },
					{ 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0 },
					{ 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0 },
					{ 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0 },
					{ 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0 },
					{ 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0 },
					{ 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0 },
					{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0 },
					{ 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0 },
					{ 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0 },
					{ 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0 },
					{ 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0 },
					{ 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 1.0 },
					{ 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0 },
					{ 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0 },
					{ 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0 },
					{ 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0 },
					{ 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0 },
					{ 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 0.0, 0.0 },
					{ 1.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0 },
					{ 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0 },
					{ 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0, 1.0 },
					{ 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0 },
			};
			String IDUM = JOptionPane.showInputDialog
					("Enter seed for RN generators use in six digits. ");
			int idum = Integer.parseInt(IDUM);
			System.out.println("RN seed used was= " + idum);
			int iseed = - idum; // test line
			Random generator2 = new Random(idum);
			double trn = generator2.nextDouble(); // test line
			double fmu = 0.0;
			double var = 1.0;
			double pma = 1.0;
			for (int j = 0; j <= NJ - 1; j++) {
				a[j] = pma;
			} // end for (j ...
			for (int j = 0; j <= NJ - 1; j++) {
				double sum1 = 0.0;
				for (int i = 0; i < NI - 1; i++) {
					sum1 += y[i][j];
				} // end for (i ...
				double phat = (sum1 + 0.5) / (NI + 1.0);
				g[j] = -phiinv(phat) * Math.sqrt(1.0 + pma);
			} // end for (j ...
			for (int i = 0; i <= NI - 1; i++) {
				th[i] = 0.0;
			} // end for (i ...
			String ito = JOptionPane.showInputDialog
					("Enter iteration number to start recording. ");
			int ITO = Integer.parseInt(ito);
			System.out.println
					("Iteration number to start recording= " + ITO);
			String iwo = JOptionPane.showInputDialog
					("Enter number of iterations between recordings. ");
			int IWO = Integer.parseInt(iwo);
			System.out.println
					("Number of iterates between recordings (t)= " + IWO);
			String num = JOptionPane.showInputDialog
					("Enter number of recorded iterates. ");
			int NUM = Integer.parseInt(num);
			System.out.println
					("Number of recorded iterates (num)= " + NUM);
			int ITOTAL = ITO + (NUM - 1) * IWO;
			int idiff = IWO - 1;
			for (int kk = 0; kk <= ITOTAL - 1; kk++) {
				System.out.println("Iteration= " + (kk + 1));
				for (int j = 0; j <= NJ - 1; j++) {
					for (int i = 0; i <= NI - 1; i++) {
						double flp = th[i] * a[j] - g[j];
						double bb = phi(-flp);
						double u = generator2.nextDouble();
						double tt = (bb * (1.0 - y[i][j]) +
								(1.0 - bb) * y[i][j]) * u + bb * y[i][j];
						z[i][j] = phiinv(tt) + flp;
					} // end for (int i ...
				} // end for (int j ...
				double suma2 = 0.0;
				for (int j = 0; j <= NJ - 1; j++) {
					suma2 += a[j] * a[j];
				} // end for (int j ...
				double v = 1.0 / suma2;
				double pvar = 1.0 / (1.0 / v + 1.0 / var);
				for (int i = 0; i <= NI - 1; i++) {
					double fmn = 0.0;
					for (int j = 0; j <= NJ - 1; j++) {
						fmn += a[j] * (z[i][j] + g[j]);
					} // end for (int j ...
					double pmean = (fmn + fmu / var) * pvar;
					th[i] = generator2.nextGaussian() * Math.sqrt(pvar) +
							pmean;
				} // end for (int i ...
				xpx[0][0] = 0.0;
				xpx[0][1] = 0.0;
				xpx[1][0] = 0.0;
				xpx[1][1] = NI;
				for (int i = 0; i <= NI - 1; i++) {
					xpx[0][0] += th[i] * th[i];
					xpx[0][1] += th[i];
				} // end for (int i ...
				xpx[0][1] = -xpx[0][1];
				xpx[1][0] = xpx[0][1];
				double sda = 0.5;
				double sdg = 2.0;
				pp[0][0] = 1.0 / Math.pow(sda, 2);
				pp[0][1] = 0.0;
				pp[1][0] = 0.0;
				pp[1][1] = 1.0 / Math.pow(sdg, 2);
				xpx[0][0] += pp[0][0];
				xpx[0][1] += pp[0][1];
				xpx[1][0] += pp[1][0];
				xpx[1][1] += pp[1][1];
				double det = xpx[0][0] * xpx[1][1] - xpx[0][1] * xpx[1][0];
				xpxinv[0][0] = xpx[1][1] / det;
				xpxinv[0][1] = -xpx[0][1] / det;
				xpxinv[1][0] = -xpx[1][0] / det;
				xpxinv[1][1] = xpx[0][0] / det;
				amat[0][0] = Math.sqrt(xpxinv[0][0]);
				amat[0][1] = xpxinv[0][1] / amat[0][0];
				amat[1][0] = 0.0;
				amat[1][1] = Math.sqrt(xpxinv[1][1] -
						(xpxinv[0][1] * xpxinv[0][1] / xpxinv[0][0]));
				for (int j = 0; j <= NJ - 1; j++) {
					double xpz1 = 0.0;
					double xpz2 = 0.0;
					for (int i = 0; i <= NI - 1; i++) {
						xpz1 += th[i] * z[i][j];
						xpz2 += z[i][j];
					} // end for (int i...
					xpz2 = -xpz2;
					pma = 1.0;
					double pmb = 0.0;
					double bz1 = xpxinv[0][0] * (xpz1 + pp[0][0] * pma) +
							xpxinv[0][1] * (xpz2 + pp[1][1] * pmb);
					double bz2 = xpxinv[1][0] * (xpz1 + pp[0][0] * pma) +
							xpxinv[1][1] * (xpz2 + pp[1][1] * pmb);
					amatt[0][0] = amat[0][0];
					amatt[0][1] = amat[1][0];
					amatt[1][0] = amat[0][1];
					amatt[1][1] = amat[1][1];
					do {
						double rn1 = generator2.nextGaussian();
						double rn2 = generator2.nextGaussian();
						a[j] = amatt[0][0] * rn1 + amatt[0][1] * rn2 + bz1;
						g[j] = amatt[1][0] * rn1 + amatt[1][1] * rn2 + bz2;
					} while (a[j] <= 0.0);
				} // end for (int j...
				if (kk >= ITO - 1) {
					idiff += 1;
					if (idiff == IWO) {
						idiff = 0;
						for (int j = 0; j <= NJ - 1; j++) {
							System.out.println("KK= " + (kk + 1) +
									"  J= " + (j + 1) +
									"  A(J)= " + a[j] +
									"  G(J)= " + g[j]);
						} // end for (int j ...
						for (int i = 0; i <= NI - 1; i++) {
							System.out.println("I= " + (i + 1) +
									"  TH(I)= " + th[i]);
						} // end for (int i ...
						System.out.println("A(J)");
						for (int j = 0; j <= NJ - 1; j++) {
							System.out.println(a[j]);
						} // end for (int j ...
						System.out.println("G(J)");
						for (int j = 0; j <= NJ - 1; j++) {
							System.out.println(g[j]);
						} // end for (int j ...
					}
					else {
						// normal processing
					} // end if (idiff ... else ...
				}
				else {
					// normal processing
				} // end if (int kk ... else ...
			} // end for (int kk ...
		} // end public static

		public static double phi(double x) {
			// (Abramowitz & Stegun, 1972, p. 932)
			double y, p;
			final double d[] = {.049867347,
					.0211410061,
					.0032776263,
					.0000380036,
					.0000488906,
					.000005383};
			if (x < -4.87726) {
				p = .00001;
				return p;
			}
			else if (x > 4.87726) {
				p = .99999;
				return p;
			}
			else {
				y = 0.;
				if (x < 0.) {
					y = -1.;
					x = -x;
				}
				else {
					// normal processing
				} // end if (x ... else ...

				p = 1. - (1. / 2. ) / Math.pow((1. +
						d[0] * x +
						d[1] * Math.pow(x, 2) +
						d[2] * Math.pow(x, 3) +
						d[3] * Math.pow(x, 4) +
						d[4] * Math.pow(x, 5) +
						d[5] * Math.pow(x, 6)) ,16);
				if (y == -1.) {
					p = (1. / 2. ) / Math.pow((1. +
							d[0] * x +
							d[1] * Math.pow(x, 2) +
							d[2] * Math.pow(x, 3) +
							d[3] * Math.pow(x, 4) +
							d[4] * Math.pow(x, 5) +
							d[5] * Math.pow(x, 6)) ,16);
				}
				else{
					// normal processing
				} // end if (y ...
				return p;
			} // end if (x ... elseif ... else ...
		} // end public static

		public static double phiinv(double p) {
			// (Abramowitz & Stegun, 1972, p. 933)
			double y, x, t;
			final double c[] = {2.515517,
					.802853,
					.010328};
			final double d[] = {1.432788,
					.189269,
					.001308};
			if (p < .00001) {
				x = -4.87762;
				return x;
			}
			else if (p > .99999) {
				x = 4.87762;
				return x;
			}
			else {
				y = 0.;
				if (p > .5) {
					y = -1.;
					p = 1. - p;
				}
				else {
					// normal processing
				} // end if (x ... else ...
				t = Math.sqrt(Math.log(1./ Math.pow(p, 2)));
				x = - t + (c[0] + c[1] * t + c[2] * Math.pow(t, 2)) /
						(1. + d[0] * t + d[1] * Math.pow(t,2) + d[2] *
								Math.pow(t,3));
				if (y == -1.) {
					x = t - (c[0] + c[1] * t + c[2] * Math.pow(t, 2)) /
							(1. + d[0] * t + d[1] * Math.pow(t,2) + d[2] *
									Math.pow(t,3));
					p = 1. - p;
				}
				else{
					// normal processing
				} // end if (y ...
				return x;
			} // end if (x ... elseif ... else ...
		} // end public static
	} // end public class ...
