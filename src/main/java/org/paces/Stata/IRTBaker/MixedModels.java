package org.paces.Stata.IRTBaker;

import javax.swing.*;

/**
 * Created by billy on 9/27/15.
 */
public class MixedModels {


/**
 * AppendixJ.java
 * Estimation of mixed models.
 */
		static final int NL = 58; // Number of theta Levels
		static final int NI = 5; // Number of items I
		static final int NQ = 10; // Number of Quadratures
		//static final int NM = 15; // No. of all paraMeters - do not use this
		static final int MK = 3; // Max(K(i)) Maximum of categories K
		static final int NS = 3; // the exact size N of the matrix inverted S
		static final int NT = 2; // the exact size N of the matrix inverted T
		static final int MNC = 25; // Maximum Number of Cycles
		static final int MIT = 10; // Maximum number of ITerations
		public MixedModels(String[] args) {
//  double[] x = new double[NQ];
			double[] ax = new double[NQ];
//  double[] r = new double[NL];
			double[] pb = new double[NL];
			double[] pa = new double[MK+1]; // n.b. subscript
			double[] p = new double[MK];
			double[] dz = new double[MK];
			double[] dzeta = new double[MK];
			double[][] a1 = new double[NS][NS];
			double[][] y1 = new double[NS][NS];
			int[] indx1 = new int[NS];
			double[][] a2 = new double[NT][NT];
			double[][] y2 = new double[NT][NT];
			int[] indx2 = new int[NT];
			// y[NL][NI]=missing indicator
//  double[][] uu = new double[NL][NI];
			double[][][] u = new double[NL][NI][MK];
//  double[][] y = new double[NL][NI];
//  double[][] zeta = new double[NI][MK-1];
//  double[] lambda = new double[NI];
			double[][] lx = new double[NL][NQ];
			double[][][] rb = new double[NI][MK][NQ];
			double[][][] ib = new double[NI][MK][NQ];
			double[][] ddzz = new double[MK][MK];
			double[][] wa = new double[MK+1][NQ]; // n.b. subscript
			double[][] px = new double[MK][NQ];
			double[] col2 = new double[MK-1];
			double[] col1 = new double[MK]; // matrix inversion
			int[] nk = {
					2, 2, 2, 3, 3
			};
			double[][] uu = {
					{ 1.0, 1.0, 1.0, 1.0, 1.0 },
					{ 1.0, 1.0, 1.0, 1.0, 3.0 },
					{ 1.0, 1.0, 1.0, 2.0, 1.0 },
					{ 1.0, 1.0, 1.0, 2.0, 2.0 },
					{ 1.0, 1.0, 1.0, 3.0, 1.0 },
					{ 1.0, 1.0, 1.0, 3.0, 2.0 },
					{ 1.0, 1.0, 2.0, 1.0, 2.0 },
					{ 1.0, 1.0, 2.0, 2.0, 1.0 },
					{ 1.0, 1.0, 2.0, 2.0, 2.0 },
					{ 1.0, 1.0, 2.0, 2.0, 3.0 },
					{ 1.0, 1.0, 2.0, 3.0, 2.0 },
					{ 1.0, 1.0, 2.0, 3.0, 3.0 },
					{ 1.0, 2.0, 1.0, 1.0, 1.0 },
					{ 1.0, 2.0, 1.0, 1.0, 2.0 },
					{ 1.0, 2.0, 1.0, 2.0, 1.0 },
					{ 1.0, 2.0, 1.0, 2.0, 2.0 },
					{ 1.0, 2.0, 1.0, 3.0, 1.0 },
					{ 1.0, 2.0, 1.0, 3.0, 2.0 },
					{ 1.0, 2.0, 1.0, 3.0, 3.0 },
					{ 1.0, 2.0, 2.0, 1.0, 1.0 },
					{ 1.0, 2.0, 2.0, 1.0, 2.0 },
					{ 1.0, 2.0, 2.0, 2.0, 1.0 },
					{ 1.0, 2.0, 2.0, 2.0, 2.0 },
					{ 1.0, 2.0, 2.0, 2.0, 3.0 },
					{ 1.0, 2.0, 2.0, 3.0, 1.0 },
					{ 1.0, 2.0, 2.0, 3.0, 2.0 },
					{ 1.0, 2.0, 2.0, 3.0, 3.0 },
					{ 2.0, 1.0, 1.0, 1.0, 1.0 },
					{ 2.0, 1.0, 1.0, 2.0, 1.0 },
					{ 2.0, 1.0, 1.0, 2.0, 2.0 },
					{ 2.0, 1.0, 1.0, 2.0, 3.0 },
					{ 2.0, 1.0, 1.0, 3.0, 1.0 },
					{ 2.0, 1.0, 1.0, 3.0, 3.0 },
					{ 2.0, 1.0, 2.0, 1.0, 1.0 },
					{ 2.0, 1.0, 2.0, 1.0, 2.0 },
					{ 2.0, 1.0, 2.0, 1.0, 3.0 },
					{ 2.0, 1.0, 2.0, 2.0, 1.0 },
					{ 2.0, 1.0, 2.0, 2.0, 2.0 },
					{ 2.0, 1.0, 2.0, 2.0, 3.0 },
					{ 2.0, 1.0, 2.0, 3.0, 2.0 },
					{ 2.0, 1.0, 2.0, 3.0, 3.0 },
					{ 2.0, 2.0, 1.0, 1.0, 1.0 },
					{ 2.0, 2.0, 1.0, 1.0, 2.0 },
					{ 2.0, 2.0, 1.0, 2.0, 1.0 },
					{ 2.0, 2.0, 1.0, 2.0, 2.0 },
					{ 2.0, 2.0, 1.0, 2.0, 3.0 },
					{ 2.0, 2.0, 1.0, 3.0, 1.0 },
					{ 2.0, 2.0, 1.0, 3.0, 2.0 },
					{ 2.0, 2.0, 1.0, 3.0, 3.0 },
					{ 2.0, 2.0, 2.0, 1.0, 1.0 },
					{ 2.0, 2.0, 2.0, 1.0, 2.0 },
					{ 2.0, 2.0, 2.0, 1.0, 3.0 },
					{ 2.0, 2.0, 2.0, 2.0, 1.0 },
					{ 2.0, 2.0, 2.0, 2.0, 2.0 },
					{ 2.0, 2.0, 2.0, 2.0, 3.0 },
					{ 2.0, 2.0, 2.0, 3.0, 1.0 },
					{ 2.0, 2.0, 2.0, 3.0, 2.0 },
					{ 2.0, 2.0, 2.0, 3.0, 3.0 }
			};
			double[][] y = {
					{ 1.0, 1.0, 1.0, 1.0, 1.0 },
					{ 1.0, 1.0, 1.0, 1.0, 1.0 },
					{ 1.0, 1.0, 1.0, 1.0, 1.0 },
					{ 1.0, 1.0, 1.0, 1.0, 1.0 },
					{ 1.0, 1.0, 1.0, 1.0, 1.0 },
					{ 1.0, 1.0, 1.0, 1.0, 1.0 },
					{ 1.0, 1.0, 1.0, 1.0, 1.0 },
					{ 1.0, 1.0, 1.0, 1.0, 1.0 },
					{ 1.0, 1.0, 1.0, 1.0, 1.0 },
					{ 1.0, 1.0, 1.0, 1.0, 1.0 },
					{ 1.0, 1.0, 1.0, 1.0, 1.0 },
					{ 1.0, 1.0, 1.0, 1.0, 1.0 },
					{ 1.0, 1.0, 1.0, 1.0, 1.0 },
					{ 1.0, 1.0, 1.0, 1.0, 1.0 },
					{ 1.0, 1.0, 1.0, 1.0, 1.0 },
					{ 1.0, 1.0, 1.0, 1.0, 1.0 },
					{ 1.0, 1.0, 1.0, 1.0, 1.0 },
					{ 1.0, 1.0, 1.0, 1.0, 1.0 },
					{ 1.0, 1.0, 1.0, 1.0, 1.0 },
					{ 1.0, 1.0, 1.0, 1.0, 1.0 },
					{ 1.0, 1.0, 1.0, 1.0, 1.0 },
					{ 1.0, 1.0, 1.0, 1.0, 1.0 },
					{ 1.0, 1.0, 1.0, 1.0, 1.0 },
					{ 1.0, 1.0, 1.0, 1.0, 1.0 },
					{ 1.0, 1.0, 1.0, 1.0, 1.0 },
					{ 1.0, 1.0, 1.0, 1.0, 1.0 },
					{ 1.0, 1.0, 1.0, 1.0, 1.0 },
					{ 1.0, 1.0, 1.0, 1.0, 1.0 },
					{ 1.0, 1.0, 1.0, 1.0, 1.0 },
					{ 1.0, 1.0, 1.0, 1.0, 1.0 },
					{ 1.0, 1.0, 1.0, 1.0, 1.0 },
					{ 1.0, 1.0, 1.0, 1.0, 1.0 },
					{ 1.0, 1.0, 1.0, 1.0, 1.0 },
					{ 1.0, 1.0, 1.0, 1.0, 1.0 },
					{ 1.0, 1.0, 1.0, 1.0, 1.0 },
					{ 1.0, 1.0, 1.0, 1.0, 1.0 },
					{ 1.0, 1.0, 1.0, 1.0, 1.0 },
					{ 1.0, 1.0, 1.0, 1.0, 1.0 },
					{ 1.0, 1.0, 1.0, 1.0, 1.0 },
					{ 1.0, 1.0, 1.0, 1.0, 1.0 },
					{ 1.0, 1.0, 1.0, 1.0, 1.0 },
					{ 1.0, 1.0, 1.0, 1.0, 1.0 },
					{ 1.0, 1.0, 1.0, 1.0, 1.0 },
					{ 1.0, 1.0, 1.0, 1.0, 1.0 },
					{ 1.0, 1.0, 1.0, 1.0, 1.0 },
					{ 1.0, 1.0, 1.0, 1.0, 1.0 },
					{ 1.0, 1.0, 1.0, 1.0, 1.0 },
					{ 1.0, 1.0, 1.0, 1.0, 1.0 },
					{ 1.0, 1.0, 1.0, 1.0, 1.0 },
					{ 1.0, 1.0, 1.0, 1.0, 1.0 },
					{ 1.0, 1.0, 1.0, 1.0, 1.0 },
					{ 1.0, 1.0, 1.0, 1.0, 1.0 },
					{ 1.0, 1.0, 1.0, 1.0, 1.0 },
					{ 1.0, 1.0, 1.0, 1.0, 1.0 },
					{ 1.0, 1.0, 1.0, 1.0, 1.0 },
					{ 1.0, 1.0, 1.0, 1.0, 1.0 },
					{ 1.0, 1.0, 1.0, 1.0, 1.0 },
					{ 1.0, 1.0, 1.0, 1.0, 1.0 }
			};
			double[] r = {
					1.0, 1.0, 2.0, 4.0, 1.0, 1.0, 3.0, 1.0, 2.0, 2.0,
					1.0, 3.0, 2.0, 3.0, 8.0, 4.0, 2.0, 6.0, 4.0, 5.0,
					1.0, 10.0, 7.0, 2.0, 7.0, 10.0, 8.0, 1.0, 1.0, 3.0,
					1.0, 4.0, 4.0, 4.0, 4.0, 2.0, 3.0, 11.0, 2.0, 8.0,
					6.0, 3.0, 13.0, 11.0, 26.0, 13.0, 12.0, 31.0, 25.0, 6.0,
					8.0, 8.0, 17.0, 60.0, 54.0, 21.0, 123.0, 179.0
			};
			// initial values
			double[][] zeta = new double[NI][MK-1];
			zeta[0][0] = 0.0;
			zeta[1][0] = 0.0;
			zeta[2][0] = 0.0;
			zeta[3][0] = 0.693;
			zeta[4][0] = 0.693;
			zeta[0][1] = 0.0;
			zeta[1][1] = 0.0;
			zeta[2][1] = 0.0;
			zeta[3][1] = -0.693;
			zeta[4][1] = -0.693;
			double[] lambda = new double[NI];
			lambda[0] = 1.0;
			lambda[1] = 1.0;
			lambda[2] = 1.0;
			lambda[3] = 1.0;
			lambda[4] = 1.0;
			double[] x = {
					-4.5, -3.5, -2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 3.5, 4.5
			};
			double mu = 0.0;
			double s2 = 1.0;
			double sax = 0.0;
			for (int l = 0; l <= NL - 1; l++) {
				for (int i = 0; i <= NI - 1; i++) {
					for (int k = 0; k <= MK - 1; k++) {
						u[l][i][k] = 0.0;
					} // end for (int k ...
				} // end for (int i ...
			} // end for (int l ...
			for (int l = 0; l <= NL - 1; l++) {
				for (int i = 0; i <= NI - 1; i++) {
					for (int k = 0; k <= MK - 1; k++) {
						if (uu[l][i] == (double) k + 1) {
							u[l][i][k] = 1.0;
						} // end if (uu[l][i] ...
//         System.out.println("L= " + (l + 1) +
//           "  I= " + (i + 1) +
//           "  K= " + (k + 1) +
//           "  U(L,I,K)= " + u[l][i][k]);
					} // end for (int k ...
				} // end for (int i ...
			} // end for (int l ...
			for (int q = 0; q <= NQ - 1; q++) {
				ax[q] = (1.0 / Math.sqrt(2.0 * Math.PI)) *
						Math.exp(-(1.0 / 2.0) *
								Math.pow(x[q] - mu, 2) / s2);
				sax += ax[q];
			} // end for (int q ...
			for (int q = 0; q <= NQ - 1; q++) {
				ax[q] /= sax;
			} // end for (int q ...
			int trall = -1; // set traces
			String yn = JOptionPane.showInputDialog
					("Trace all? Y/N");
			if (yn.equalsIgnoreCase("y")) {
				trall = 0;
			}
			else{
				// normal processing
			} // end if (yn ... else
			if (trall == 0) {
				for (int q = 0; q <= NQ - 1; q++) {
					System.out.println("Q= " + (q + 1) +
							"  X= " + x[q] +
							"  AX= " + ax[q]);
				} // end for (int q ...
			}
			else {
				// normal processing
			} // end if (trall ... else
			for (int nc = 0; nc <= MNC - 1; nc++) {
				System.out.println("Cycle= " + (nc + 1));
				for (int l = 0; l <= NL - 1; l++) {
					pb[l] = 0.0;
					for (int q = 0; q <= NQ - 1; q++) {
						lx[l][q] = 1.0;
					} // end for (int q ...
				} // end for (int l ...
				// e step
				for (int l = 0; l <= NL - 1; l++) {
					for (int i = 0; i <= NI - 1; i++) {
						for (int q = 0; q <= NQ - 1; q++) {
							if (y[l][i] != 0.0) {
								// y not employed
								pa[0] = 1.0;
								pa[nk[i]] = 0.0;
								// very careful about the subscript
								for (int k = 0; k <= nk[i] - 2; k++) {
									pa[k+1] = 1.0 /
											(1.0 + Math.exp(-(zeta[i][k] + lambda[i] * x[q])));
								} // end for (int k ...
								double pq = 0.0;
								for (int k = 0; k <= nk[i] - 1; k++) {
									p[k] = pa[k] - pa[k+1];
									pq += u[l][i][k]*p[k];
									if (trall == 0) {
										if (l == 0) {
											System.out.println("I= " + (i + 1) +
													"  Q= " + (q + 1) +
													"  P(K)= " + p[k] +
													"  PQ= " + pq);
										}
										else {
											// normal processing
										} // end if (l ...
									}
									else {
										// normal processing
									} // end if (trall ... else
								} // end for (int k ...
								lx[l][q] *= pq;
							}
							else {
								// normal processing
							} // end if (y[l][i] ... else
						} // end for (int q ...
					} // end for (int i ...
					pb[l] = 0.0;
					for (int q = 0; q <= NQ - 1; q++) {
						pb[l] += lx[l][q] * ax[q];
					} // end for (int q ...
				} // end for (int l ...
				if (trall == 0) {
					for (int l = 0; l <= NL - 1; l++) {
						System.out.println("L= " + (l + 1) +
								"  PB(L)= " + pb[l]);
					} // end for (int l ...
				}
				else {
					// normal processing
				} // end if (trall ... else
				if (trall == 0) {
					for (int ll = 0; ll <= NL - 1; ll++) {
						for (int q = 0; q <= NQ - 1; q++) {
							System.out.println("L= " + (ll + 1) +
									"  Q= " + (q + 1) +
									"  LX(L,Q)= " + lx[ll][q]);
						} // end for (int q ...
					} // end for (int l ...
				}
				else {
					// normal processing
				} // end if (trall ... else
				// n bar and r bar
				for (int i = 0; i <= NI - 1; i++) {
					for (int k = 0; k <= nk[i] - 1; k++) {
						for (int q = 0; q <= NQ - 1; q++) {
							rb[i][k][q] = 0.0;
							ib[i][k][q] = 0.0;
							for (int ll = 0; ll <= NL - 1; ll++) {
								if (y[ll][i] != 0.0) {
									rb[i][k][q] += u[ll][i][k] * r[ll] * lx[ll][q]
											* ax[q] / pb[ll];
									ib[i][k][q] += r[ll] * lx[ll][q] * ax[q] /
											pb[ll];
								}
								else {
									// normal processing
								} // end if (y[g][l][i] ... else
							} // end for (int l ...
						} // end for (int q ...
					} // end for (int k ...
				} // end for (int i ...
				if (trall == 0) {
					for (int i = 0; i <= NI - 1; i++) {
						for (int k = 0; k <= nk[i] - 1; k++) {
							for (int q = 0; q <= NQ - 1; q++) {
								System.out.println("I= " + (i + 1) +
										"  K= " + (k + 1) +
										"  Q= " + (q + 1) +
										"  RB(I,K,Q)= " + rb[i][k][q] +
										"  IB(I,K,Q)= " + ib[i][k][q]);
							} // end for (int q ...
						} // end for (int k ...
					} // end for (int i ...
				}
				else {
					// normal processing
				} // end if (trall ... else
				// probit stage m step
				for (int i = 0; i <= NI - 1; i++) {
					for (int it = 0; it <= MIT - 1; it++) {
						for (int k = 0; k <= nk[i] - 1; k++) {
							dz[k] = 0.0;
							for (int kk = 0; kk <= nk[i] - 1; kk++) {
								ddzz[k][kk] = 0.0;
							} // end for (int kk ...
						} // end for (int k ...
						for (int k = 0; k <= nk[i] - 1; k++) {
							for (int q = 0; q <= NQ - 1; q ++) {
								pa[0] = 1.0;
								pa[nk[i]] = 0.0;
								wa[0][q] = 0.0;
								wa[nk[i]][q] = 0.0;
								for (int kk = 0; kk <= nk[i] - 2; kk++) {
									pa[kk+1] = 1.0 / (1.0 +
											Math.exp(-(zeta[i][kk] + lambda[i]*x[q])));
								} // end for (int kk ...
								px[k][q] = pa[k] - pa[k+1];
								wa[k+1][q] = pa[k+1]*(1.0 - pa[k+1]);
							} // end for (int q ...
						} // end for (int k ...
						for (int k = 0; k <= nk[i] - 2; k++) {
							for (int q = 0; q <= NQ - 1; q++) {
								dz[k] += (-rb[i][k][q] / px[k][q] +
										rb[i][k+1][q] / px[k+1][q])*wa[k+1][q];
							} // end for (int q ...
						} // end for (int k ...
						for (int q = 0; q <= NQ - 1; q++) {
							for (int kk = 0; kk <= nk[i] - 1; kk++) {
								dz[nk[i]-1] += (rb[i][kk][q] / px[kk][q]) *
										(wa[kk][q] - wa[kk+1][q])*x[q];
							} // end for (int kk ...
						} // end for (int q ...
						// BASIC line 445
						for (int k = 0; k <= nk[i] - 2; k++) {
							for (int q = 0; q <= NQ - 1; q++) {
								ddzz[k][k] -= ib[i][k][q] * wa[k+1][q] * wa[k+1][q] *
										(1.0 / px[k][q] + 1.0 / px[k+1][q]);
							} // end for (int k ...
						} // end for (int k ...
						for (int q = 0; q <= NQ - 1; q++) {
							for (int kk = 0; kk <= nk[i] - 1; kk++) {
								ddzz[nk[i]-1][nk[i]-1] -= ib[i][kk][q] *
										x[q] * x[q] * (wa[kk][q] - wa[kk+1][q]) *
										(wa[kk][q] - wa[kk+1][q]) / px[kk][q];
							} // end for (int kk ...
						} // end for (int q ...
						if (nk[i] >= 3) {
							for (int k = 1; k <= nk[i] - 2; k++) {
								for (int q = 0; q <= NQ - 1; q++) {
									ddzz[k][k-1] += ib[i][k][q] *
											wa[k+1][q] * wa[k][q] / px[k][q];
								} // end for (int q ...
							} // end for (int k ...
						} // end if (nk[i] ...
						if (nk[i] >= 3) {
							for (int k = 0; k <= nk[i] - 3; k++) {
								for (int q = 0; q <= NQ - 1; q++) {
									ddzz[k][k+1] += ib[i][k][q] *
											wa[k+1][q] * wa[k+2][q] / px[k+1][q];
								} // end for (int q ...
							} // end for (int k ...
						} // end if (nk[i] ...
						for (int k = 0; k <= nk[i] - 2; k++) {
							for (int q = 0; q <= NQ - 1; q++) {
								ddzz[k][nk[i]-1] += ib[i][k][q] * wa[k+1][q] * x[q] *
										((wa[k][q] - wa[k+1][q]) / px[k][q] -
												(wa[k+1][q] - wa[k+2][q]) / px[k+1][q]);
							} // end for (int q ...
							ddzz[nk[i]-1][k] = ddzz[k][nk[i]-1];
						} // end for (int k ...
						if (nk[i] == 2) {
							for (int k = 0; k <= nk[i] - 1; k++) {
								for (int kk = 0; kk <= nk[i] - 1; kk++) {
									a2[k][kk] = ddzz[k][kk];
									if (trall == 0) {
										System.out.println("K= " + (k + 1) +
												"  KK= " + (kk + 1) +
												"  A2(K,KK)= " + a2[k][kk]);
									} // if (trall == ...
								} // end for (int kk ...
							} // end for (int k ...
						}
						else {
							for (int k = 0; k <= nk[i] - 1; k++) {
								for (int kk = 0; kk <= nk[i] - 1; kk++) {
									a1[k][kk] = ddzz[k][kk];
									if (trall == 0) {
										System.out.println("K= " + (k + 1) +
												"  KK= " + (kk + 1) +
												"  A2(K,KK)= " + a1[k][kk]);
									} // end if (trall ...
								} // end for (int kk ...
							} // end for (int k ...
						} // end if (nk[i] ... else ...
						int np = nk[i]; // np not really used
						int N = np;
						if (nk[i] == 2) {
							ludcmp(a2, N, indx2);
							for (int ii = 0; ii <= N - 1; ii++) {
								for (int ll = 0; ll <= N - 1; ll++ ) col2[ll] = 0.0;
								col2[ii] = 1.0;
								lubksb(a2, N, indx2, col2);
								for (int ll = 0; ll <= N - 1; ll++) y2[ll][ii] = col2[ll];
							} // end for (int ii ...
						}
						else {
							ludcmp(a1, N, indx1);
							for (int ii = 0; ii <= N - 1; ii++) {
								for (int ll = 0; ll <= N - 1; ll++ ) col1[ll] = 0.0;
								col1[ii] = 1.0;
								lubksb(a1, N, indx1, col1);
								for (int ll = 0; ll <= N - 1; ll++) y1[ll][ii] = col1[ll];
							} // end for (int ii ...
						} // end if (nk[i] ... else ...
						for (int k = 0; k <= nk[i] - 1; k++) {
							dzeta[k] = 0.0;
						} // end for (int k ...
						if (nk[i] == 2) {
							for (int k = 0; k <= nk[i] - 1; k++) {
								for (int kk = 0; kk <= nk[i] - 1; kk++) {
									dzeta[k] += y2[k][kk] * dz[kk];
								} // end for (int kk ...
							} // end for (int k ...
						}
						else {
							for (int k = 0; k <= nk[i] - 1; k++) {
								for (int kk = 0; kk <= nk[i] - 1; kk++) {
									dzeta[k] += y1[k][kk] * dz[kk];
								} // end for (int kk ...
							} // end for (int k ...
						} // end if (nk[k] ...
						// big ezeta[k] protection
						for (int k = 0; k <= nk[i] - 1; k++) {
							for (int kk = 0; kk <= nk[i] - 1; kk++) {
								if (dzeta[k] > 0.2) dzeta[k] = 0.2;
								if (dzeta[k] < -0.2) dzeta[k] = -0.2;
							} // end for (int kk ...
						} // end for (int k ...
						for (int k = 0; k <= nk[i] - 2; k++) {
							zeta[i][k] -= dzeta[k];
						} // end for (int k ...
						lambda[i] -= dzeta[nk[i] - 1];
						System.out.println("I= " + (i + 1) +
								"  IT= " + (it + 1) +
								"  ZETA(I,1)= " + zeta[i][0] +
								"  ZETA(I,2)= " + zeta[i][1] +
								"  LAMBDA(I)= " + lambda[i]);
						double diffzl = 0.0;
						for (int k = 0; k <= nk[i] - 1; k ++ ) {
							if (Math.abs(dzeta[k]) > diffzl) {
								diffzl = Math.abs(dzeta[k]);
							}
							else {
								// normal processing
							} // end if (Math.abs(dzeta[k]) ...
						} // end for (int k ...
						if (diffzl <= .01) {
							break; // exit for (int it ...
						}
						else {
							// normal processing
						} // end if (diffzl ...
					} // end for (int it ... BASIC line 420
				} // end for (int i ...
				for (int i = 0; i <= NI - 1; i ++) {
					if (nk[i] == 2) {
						System.out.println("I= " + (i + 1) +
								"  ZETA(I,1)= " + zeta[i][0] +
								"  LAMBDA(I)= " + lambda[i] +
								"  B(I,1)= " + (-zeta[i][0]/lambda[i]));
					} // end if (nk[i] ...
					if (nk[i] == 3) {
						System.out.println("I= " + (i + 1) +
								"  ZETA(I,1)= " + zeta[i][0] +
								"  ZETA(I,2)= " + zeta[i][1] +
								"  LAMBDA(I)= " + lambda[i] +
								"  B(I,1)= " + (-zeta[i][0]/lambda[i]) +
								"  B(I,2)= " + (-zeta[i][1]/lambda[i]));
					} // end if (nk[i] ...
				} // end for (int i ... BASIC line 400
				// this should be modified.
			} // end for (int nc ...
		} // end public static ...
		// LU decomposition.
		public static double ludcmp(double[][] a, int n, int[] index) {
			int i, imax = 0, j, k;
			double d, big, dum, sum, temp;
			double[] vv = new double[n];

			d = 1.0;
			for (i = 0; i < n; i++) {
				big = 0.0;
				for (j = 0; j < n; j++)
					if ((temp = Math.abs(a[i][j])) > big) big = temp;
				if (big == 0.0) {
					System.out.println("Singular matrix in routine ludcmp");
					System.exit(1);
				}
				vv[i] = 1.0 / big;
			}

			for (j = 0; j < n; j++) {
				for (i = 0; i < j; i++) {
					sum = a[i][j];
					for (k = 0; k < i; k++) sum -= a[i][k] * a[k][j];
					a[i][j] = sum;
				}
				big = 0.0;
				for (i = j; i < n; i++) {
					sum = a[i][j];
					for (k = 0; k < j; k++) sum -= a[i][k] * a[k][j];
					a[i][j] = sum;
					if ((dum = vv[i] * Math.abs(sum)) >= big) {
						big = dum;
						imax = i;
					}
				}
				if (j != imax) {
					for (k = 0; k < n; k++) {
						dum = a[imax][k];
						a[imax][k] = a[j][k];
						a[j][k] = dum;
					}
					d = -d;
					vv[imax] = vv[j];
				}
				index[j] = imax;
				if (a[j][j] == 0.0) a[j][j] = 1.0e-20;
				if (j != n - 1) {
					dum = 1.0 / a[j][j];
					for (i = j + 1; i < n; i++) a[i][j] *= dum;
				}
			}

			return d;
		}
		// LU backward substitution.
		public static void lubksb(double[][] a, int n, int[] index, double b[]) {
			int i, ii = -1, ip, j;
			double sum;

			for (i = 0; i < n; i++) {
				ip = index[i];
				sum = b[ip];
				b[ip] = b[i];
				if (ii != -1) for (j = ii; j <= i - 1; j++) sum -= a[i][j] * b[j];
				else if (sum != 0.0) ii = i;
				b[i] = sum;
			}

			for (i = n - 1; i >= 0; i--) {
				sum = b[i];
				for (j = i + 1; j < n; j++) sum -= a[i][j] * b[j];
				b[i] = sum / a[i][i];
			}
		}
	} // end public class ...
