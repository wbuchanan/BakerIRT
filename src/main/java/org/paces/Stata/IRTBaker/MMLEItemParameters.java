package org.paces.Stata.IRTBaker;

import javax.swing.*;

/**
 * Created by billy on 9/27/15.
 */
public class MMLEItemParameters {

/**
 * A program to implement MMLE/EM for LSAT-6 data set.
 */
		static final int NITEM = 5;
		static final int NP = 32; // Number of Patterns
		static final int NQ = 10; // Number of Quadratures
		static final int NI = 6; // Number of Iterations
		public MMLEItemParameters(String[] args) {
			double[][] u = { // read canned data
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
			// quadrature points and weights via BILOG
			double[] x = { -4.000, -3.111, -2.222, -1.333, -0.4444,
					0.4444, 1.333, 2.222, 3.111, 4.000 };
			double[] ak = { 0.000119, 0.002805, 0.03002, 0.1458, 0.3213,
					0.3213, 0.1458, 0.03002, 0.002805, 0.000119 };
			double[] fpt = { 3.0, 6.0, 2.0, 11.0, 1.0, 1.0, 3.0, 4.0, 1.0,
					8.0, 0.0, 16.0, 0.0, 3.0, 2.0, 15.0, 10.0, 29.0, 14.0, 81.0,
					3.0, 28.0, 15.0, 80.0, 16.0, 56.0, 21.0, 173.0, 11.0, 61.0,
					28.0, 298.0 };
			double[] cpt = { 0.0, 0.0, 0.0, 0.0, 0.0 };
			double[] a = { 1.0, 1.0, 1.0, 1.0, 1.0 };
			double[][] lxk = new double[NP][NQ];
			double[] pl = new double[NP];
			double[][] rik = new double[NITEM][NQ];
			double[][] nk = new double[NITEM][NQ];
			double[] pq = new double[NQ];
			double[] b = new double[NITEM];
			int trall = -1; // set up traceing
			int trlp = -1;
			int trrn = -1;
			int trmle = -1;
			String yn = JOptionPane.showInputDialog
					("Trace all? Y/N");
			if (yn.equalsIgnoreCase("y")) {
				trall = 0;
				trlp = 0;
				trrn = 0;
				trmle = 0;
			}
			else{
				String ynlp = JOptionPane.showInputDialog
						("Trace L,P? Y/N");
				if (ynlp.equalsIgnoreCase("y")) {
					trlp = 0;
				}
				String ynrn = JOptionPane.showInputDialog
						("Trace R and N? Y/N");
				if (ynrn.equalsIgnoreCase("y")) {
					trrn = 0;
				}
				String ynmle = JOptionPane.showInputDialog
						("Trace MLE? Y/N");
				if (ynmle.equalsIgnoreCase("y")) {
					trmle = 0;
				}
			} // end if (yn ... else
			String smnc = JOptionPane.showInputDialog
					("Enter number of cycles to do. ");
			int mnc = Integer.parseInt(smnc);
			System.out.println("All= " + trall +
					"  LP= " + trlp +
					"  TRRN= " + trrn +
					"  MLE= " + trmle +
					"  Max Cycles= " + mnc);
			for (int nc = 1; nc <= mnc; nc++) {
				System.out.println("Cycle= " + nc);
				for (int l = 0; l <= NP - 1; l++) {
					pl[l] = 0.0;
					for (int k = 0; k <= NQ - 1; k++) {
						lxk[l][k] = 1.0;
					} // end for (int l ...
				} // end for (int l ...
				// E-step get exptected R and N
				for (int l = 0; l <= NP - 1; l++) { // E-step subroutine
					for (int k = 0; k <= NQ - 1; k++) {
						for (int i = 0; i <= NITEM - 1; i++) {
							// subroutine to compute p(x)
							double dev = - (cpt[i] + a[i] * x[k]);
							double ep = Math.exp(dev);
							double p = 1.0 / (1.0 + ep);
							double q = 1.0 - p;
							pq[i] = p;
							if (u[l][i] == 0.0) {
								pq[i] = q;
							}
							if (trall == 0) {
								System.out.println("345 L= " + (l + 1) +
										"  K= " + (k + 1) +
										"  I= " + (i + 1) +
										"  PQ(I)= " + pq[i] +
										"  U(L,I)= " + u[l][i]);
							}
							lxk[l][k] *= pq[i];
						} // end for (int i ... (BASIC 370)
						if (trall == 0 || trlp == 0) {
							System.out.println("L= " + (l + 1) +
									"  K= " + (k + 1) +
									"  L(X)= " + lxk[l][k]);
						}
					} // end for (int k ... (BASIC 380)
					pl[l] = 0.0;
					for (int k = 0; k <= NQ - 1; k++) {
						pl[l] += lxk[l][k] * ak[k];
						if (trall == 0) {
							System.out.println("392 L= " + (l + 1) +
									"  PL(L)= " + pl[l] +
									"  LXK(L,K)= " + lxk[l][k] +
									"  AK(K)= " + ak[k]);
						}
					} // end for (int k ... (BASIC 395)
					if (trall == 0 || trlp == 0) {
						System.out.println("PL(" + (l + 1) + ")= " + pl[l]);
					}
				} // end for (int l ... (BASIC 400)
				// R-bar and N-bar loop
				for (int i = 0; i <= NITEM - 1; i++) {
					if (trall == 0 || trrn == 0) {
						System.out.println("Item= " + (i + 1));
					}
					for (int k = 0; k <= NQ - 1; k++) {
						rik[i][k] = 0.0;
						nk[i][k] = 0.0;
						for (int l = 0; l <= NP - 1; l++) {
							double nt = fpt[l] * lxk[l][k] * ak[k] / pl[l];
							double rt = nt * u[l][i];
							rik[i][k] += rt;
							nk[i][k] += nt;
							if (trall == 0) {
								System.out.println("442 I= " + (i + 1) +
										"  K= " + (k + 1) +
										"  L= " + (l + 1) +
										"  LXK(L,K)= " + lxk[l][k] +
										"  AK(K)= " + ak[k] +
										"  PL(L)= " + pl[l]);
								System.out.println("444 U(L,I)= " + u[l][i] +
										"  NT= " + nt +
										"  RT= " + rt);
								System.out.println("446 RIK(I,K)= " + rik[i][k] +
										"  NK(I,K)= " + nk[i][k]);
							}
						} // end for (int l ... (BASIC 448)
						if (trall == 0 || trrn == 0) {
							System.out.println("449 X(" + (k + 1) + ")= " + x[k] +
									"  R= " + rik[i][k] +
									"  N= " + nk[i][k]);
						}
					} // end for (int k ... (BASIC 450)
				} // end for (int i ... (BASIC 460)
				// M-step MLE of item parameters
				for (int i = 0; i <= NITEM - 1; i++) {
					// itembio routine for estimation of item parameters
					search:
					for (int nit = 1; nit <= NI; nit++) {
						double snw = 0.0;
						double snwv = 0.0;
						double snwx = 0.0;
						double snwxv = 0.0;
						double snwx2 = 0.0;
						double snwv2 = 0.0;
						for (int k = 0; k <= NQ - 1; k++) {
							if (nk[i][k] == 0.0) {
								break; // exit if (nk ...
							}
							else {
								double pi = rik[i][k] / nk[i][k];
								// calculate Newton-Raphson terms
								double dv = cpt[i] + a[i] * x[k];
								double ph = 1.0 / (1.0 + Math.exp(-dv));
								double w = ph * (1.0 - ph);
								if (w < 0.0000009) {
									break; // exit if (w ...
								}
								else {
									double v = (pi - ph) / w;
									if (trall == 0) {
										System.out.println("X(" + (k + 1) + ")= " + x[k] +
												"  PI= " + pi +
												"  P-HAT= " + ph);
										System.out.println("W= " + w +
												"  V= " + v);
									} // end if (trall ...
									double p1 = nk[i][k] * w;
									double p2 = p1 * v;
									double p3 = p2 * v;
									double p4 = p1 * x[k];
									double p5 = p4 * x[k];
									double p6 = p4 * v;
									snw += p1;
									snwv += p2;
									snwx += p4;
									snwxv += p6;
									snwx2 += p5;
									snwv2 += p3;
									if (trall == 0) {
										System.out.println("P1= " + p1 +
												"  P2= " + p2 +
												"  P3= " + p3);
										System.out.println("P4= " + p4 +
												"  P5= " + p5 +
												"  P6= " + p6);
										System.out.println("SNW= " + snw +
												"  SNWV= " + snwv +
												"  SNWX= " + snwx);
										System.out.println("SNWXV= " + snwxv +
												"  SNWX2= " + snwx2 +
												"  SNWV2= " + snwv2);
									} // end if (trall ...
								} // end if (w ...
							} // end if (nk ...
						} // end for (int k ... (BASIC 1180)
						if (snw <= 0.0) {
							System.out.println(
									"Out of bounds error: Iteration= " + nit);
							break search;
						}
						else {
							double dm = snw * snwx2 - snwx * snwx;
							if (trall == 0) {
								System.out.println("Denominator= " + dm);
							}
							if (dm <= 0.000099) {
								System.out.println(
										"Denominator error: Iteration= " + nit);
								break search;
							}
							double dcpt = (snwv * snwx2 - snwxv * snwx) / dm;
							double da = (snw * snwxv - snwx * snwv) / dm;
							cpt[i] += dcpt;
							a[i] += da;
							if (trall == 0 || trmle == 0) {
								System.out.println("NIT= " + nit +
										"  Item= " + (i + 1));
								System.out.println("Intercept= " + cpt[i] +
										"  Change= " + dcpt);
								System.out.println("Slope= " + a[i] +
										"  Change= " + da);
							}
							if (Math.abs(cpt[i]) > 30.0 ||
									Math.abs(a[i]) > 20.0) {
								System.out.println(
										"CPT(I) or A(I) error: Iteration= " + nit);
								break search;
							}
							if (Math.abs(dcpt) <= 0.05 &&
									Math.abs(da) <= 0.05) {
								break;
							}
						} // end if (snw ...
					} // end for (int nit ... (BASIC 1180)
				} // end for (int i ... (BASIC 210)
			} // end for (int nc ... (BASIC 220)
			System.out.println("Done Max Cycles");
			// print item parameter estimates
			for (int i = 0; i <= NITEM - 1; i++) {
				b[i] = - cpt[i] / a[i];
				System.out.println("Item= " + (i + 1) +
						"  Intercept= " + cpt[i] +
						"  Slope= " + a[i] +
						"  Difficulty= " + b[i]);
			} // end for (int i ... (BASIC 1740)
			System.exit(0); // last line
		} // end public static ...
	} // end public class ...
