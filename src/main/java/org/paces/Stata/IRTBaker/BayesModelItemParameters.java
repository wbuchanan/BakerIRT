package org.paces.Stata.IRTBaker;

import javax.swing.*;

/**
 * A program to implement marginalized Bayesian item parameter
 * estimation using b=difficulty, a=discrimination for LSAT-6
 * data.
 */
public class BayesModelItemParameters {

		static final int NITEM = 5;
		static final int NP = 32; // Number of Patterns
		static final int NQ = 10; // Number of Quadratures
		static final int NI = 6; // Number of Iterations
		static final double BIGT = 0.5;
		static final double HPM = 0.0; // set hyperparameters
		static final double HPVAR = 1.0;


	public BayesModelItemParameters(String[] args) {
		this.marginalizedBayesianItemParameters();
		this.modelBayesianItemParameters();
		this.bayesianAPrioriPersonParameter();
	}

		public void marginalizedBayesianItemParameters() {
			double[][] u = { // read canned data
					{ 0.0, 0.0, 0.0, 0.0, 0.0 },
					{ 0.0, 0.0, 0.0, 0.0, 1.0 },
					{ 0.0, 0.0, 0.0, 1.0, 0.0 },
					{ 0.0, 0.0, 0.0, 1.0, 1.0 },
					{ 0.0, 0.0, 1.0, 0.0, 0.0 },
					{ 0.0, 0.0, 1.0, 0.0, 1.0 },
					{ 0.0, 0.0, 1.0, 1.0, 0.0 },
					{ 0.0, 0.0, 1.0, 1.0, 1.0 },
					{ 0.0, 1.0, 0.0, 0.0, 0.0 },
					{ 0.0, 1.0, 0.0, 0.0, 1.0 },
					{ 0.0, 1.0, 0.0, 1.0, 0.0 },
					{ 0.0, 1.0, 0.0, 1.0, 1.0 },
					{ 0.0, 1.0, 1.0, 0.0, 0.0 },
					{ 0.0, 1.0, 1.0, 0.0, 1.0 },
					{ 0.0, 1.0, 1.0, 1.0, 0.0 },
					{ 0.0, 1.0, 1.0, 1.0, 1.0 },
					{ 1.0, 0.0, 0.0, 0.0, 0.0 },
					{ 1.0, 0.0, 0.0, 0.0, 1.0 },
					{ 1.0, 0.0, 0.0, 1.0, 0.0 },
					{ 1.0, 0.0, 0.0, 1.0, 1.0 },
					{ 1.0, 0.0, 1.0, 0.0, 0.0 },
					{ 1.0, 0.0, 1.0, 0.0, 1.0 },
					{ 1.0, 0.0, 1.0, 1.0, 0.0 },
					{ 1.0, 0.0, 1.0, 1.0, 1.0 },
					{ 1.0, 1.0, 0.0, 0.0, 0.0 },
					{ 1.0, 1.0, 0.0, 0.0, 1.0 },
					{ 1.0, 1.0, 0.0, 1.0, 0.0 },
					{ 1.0, 1.0, 0.0, 1.0, 1.0 },
					{ 1.0, 1.0, 1.0, 0.0, 0.0 },
					{ 1.0, 1.0, 1.0, 0.0, 1.0 },
					{ 1.0, 1.0, 1.0, 1.0, 0.0 },
					{ 1.0, 1.0, 1.0, 1.0, 1.0 }
			};
			// quadrature points and weights via BILOG
			double[] x = { -4.000, -3.111, -2.222, -1.333, -0.4444,
					0.4444, 1.333, 2.222, 3.111, 4.000 };
			double[] ak = { 0.000119, 0.002805, 0.03002, 0.1458, 0.3213,
					0.3213, 0.1458, 0.03002, 0.002805, 0.000119 };
			double[] fpt = { 3.0, 6.0, 2.0, 11.0, 1.0, 1.0, 3.0, 4.0, 1.0,
					8.0, 0.0, 16.0, 0.0, 3.0, 2.0, 15.0, 10.0, 29.0, 14.0, 81.0,
					3.0, 28.0, 15.0, 80.0, 16.0, 56.0, 21.0, 173.0, 11.0, 61.0,
					28.0, 298.0 }; // item score pattern frequencies
			double[] b = { 0.0, 0.0, 0.0, 0.0, 0.0 }; // initial values
			double[] a = { 1.0, 1.0, 1.0, 1.0, 1.0 };
			double[] pm = { 0.0, 0.0, 0.0, 0.0, 0.0 }; // prior mean
			double[] psd = { 0.5, 0.5, 0.5, 0.5, 0.5 }; // prior sd
			double[] pv = new double[NITEM]; // prior variance on log(a)
			for (int i = 0; i <= NITEM - 1; i++) {
				pv[i] = psd[i] * psd[i];
			} // end for (int i ...
			double[][] lxk = new double[NP][NQ];
			double[] pl = new double[NP];
			double[][] rik = new double[NITEM][NQ];
			double[][] nk = new double[NITEM][NQ];
			double[] pq = new double[NQ];
			double[] alpha = new double[NITEM];
			int trall = -1; // set up traceing
			int trlp = -1;
			int trrn = -1;
			int trbme = -1;
			String yn = JOptionPane.showInputDialog
					("Trace all? Y/N");
			if (yn.equalsIgnoreCase("y")) {
				trall = 0;
				trlp = 0;
				trrn = 0;
				trbme = 0;
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
				String ynbme = JOptionPane.showInputDialog
						("Trace BME? Y/N");
				if (ynbme.equalsIgnoreCase("y")) {
					trbme = 0;
				}
			} // end if (yn ... else
			int foflg = -1;
			String ynflo = JOptionPane.showInputDialog
					("Use FLOAT option? Y/N");
			if (ynflo.equalsIgnoreCase("y")) {
				foflg = 0;
			}
			System.out.println("FLOAT option= " + foflg);
			String smnc = JOptionPane.showInputDialog
					("Enter number of cycles to do. ");
			int mnc = Integer.parseInt(smnc);
			System.out.println("All= " + trall +
					"  LP= " + trlp +
					"  TRRN= " + trrn +
					"  BME= " + trbme +
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
							double dev = - a[i] * (x[k] - b[i]);
							double p = 1.0 / (1.0 + Math.exp(dev));
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
				// M-step Bayes modal estimation of item parameters
				for (int i = 0; i <= NITEM - 1; i++) {
					// Bayes modal item parameter estimation
					search:
					for (int nit = 1; nit <= NI; nit++) {
						alpha[i] = Math.log(a[i]);
						double l1 = 0.0;
						double l2 = 0.0;
						double l11 = 0.0;
						double l22 = 0.0;
						double l12 = 0.0;
						for (int k = 0; k <= NQ - 1; k++) { // theta loop
							if (nk[i][k] == 0.0) {
								break; // exit if (nk ...
							}
							else {
								if (trall == 0) {
									System.out.println("M-step Iteration= " + nit +
											"  Item= " + (i + 1) +
											"  Q-Node= " + (k + 1));
								}
								// calc P and Q
								double dv = -a[i] * (x[k] - b[i]);
								double ph = 1.0 / (1.0 + Math.exp(dv)); // p hat
								double qh = 1.0 - ph; // q hat
								double w = ph * qh;
								if (w < 0.0000009) {
									break; // exit if (w ...
								}
								else {
									double rmn = rik[i][k] - nk[i][k] * ph;
									double xmb = x[k] - b[i];
									if (trall == 0) {
										System.out.println("P= " + ph +
												"  W= " + w);
										System.out.println("R-N*P= " + rmn +
												"  X-B= " + xmb);
									} // end if (trall ...
									double p1 = nk[i][k] * w;
									double p2 = p1 * xmb;
									double p3 = p2 * xmb;
									double p4 = rmn * xmb;
									l1 += p4;
									l2 += rmn;
									l11 += p3;
									l22 += p1;
									l12 += p2;
									if (trall == 0) {
										System.out.println("P1= " + p1 +
												"  P2= " + p2);
										System.out.println("P3= " + p3 +
												"  P4= " + p4);
										System.out.println("L1= " + l1 +
												"  L2= " + l2);
										System.out.println("L11= " + l11 +
												"  L22= " + l22 +
												"  L12= " + l12);
									} // end if (trall ...
								} // end if (w ...
							} // end if (nk ...
						} // end for (int k ... (BASIC 1180)
						double ea = Math.exp(alpha[i]);
						double ea2 = ea * ea;
						double labar = pm[i];
						if (foflg == 0) {
							double sumla = 0.0;
							for (int ij = 0; ij <= NITEM - 1; ij++) {
								sumla += alpha[ij];
							} // end for (int i ...
							labar = sumla / NITEM; // casting
						}
						else {
							// normal processing
						} // end if (foflg ...
						// multiply by ea and add prior terms
						l1 = ea * l1 - (alpha[i] - labar) / pv[i];
						l2 = - ea * l2; // BASIC 1195
						l11 = ea2 * l11 + 1 / pv[i];
						l22 = ea2 * l22;
						l12 = - ea2 * l12; // BASIC 1196
						// - 1 has been introduced where appropriate in the
						// expressions of lines 1195 and 1196
						double dm = l11 * l22 - l12 * l12;
						if (trall == 0) {
							System.out.println("EA= " + ea +
									"  EA2= " + ea2);
							System.out.println("After Multiplying by " +
									"exp(alpha) and adding Bayesian terms");
							System.out.println("L1= " + l1 +
									"  L2= " + l2);
							System.out.println("L11= " + l11 +
									"  L22= " + l22 +
									"  L12= " + l12);
							System.out.println("Denominator= " + dm);
						} // end if (trall ...
						if (dm <= 0.000099) {
							System.out.println(
									"Out of bounds error: Iteration= " + nit);
							a[i] = Math.exp(alpha[i]);
							break search;
						}
						double da = (l1 * l22 - l2 * l12) / dm;
						double db = (- l1 * l12 + l2 * l11) / dm;
						alpha[i] += da;
						b[i] += db;
						if (trall == 0 || trbme == 0) {
							System.out.println("NIT= " + nit +
									"  Item= " + (i + 1));
							System.out.println("ALPHA(I)= " + alpha[i] +
									"  DA= " + da);
							System.out.println("B(I)= " + b[i] +
									"  DB= " + db);
						}
						if (Math.abs(alpha[i]) > 30.0 ||
								Math.abs(b[i]) > 20.0) {
							System.out.println(
									"ALPHA(I) or B(I) error: Iteration= " + nit);
							break search;
						}
						if (Math.abs(da) <= 0.05 &&
								Math.abs(db) <= 0.05) {
							break;
						}
					} // end for (int nit ... (BASIC 1180)
					a[i] = Math.exp(alpha[i]);
				} // end for (int i ... (BASIC 210)
			} // end for (int nc ... (BASIC 220)
			System.out.println("Done Max Cycles");
			// print item parameter estimates
			for (int i = 0; i <= NITEM - 1; i++) {
				System.out.println("Item= " + (i + 1) +
						"  Difficulty= " + b[i] +
						"  Discrimination= " + a[i]);
			} // end for (int i ... (BASIC 280)
			System.exit(0); // last line
		} // end public static ...


/**
 * A program for Bayesian modal estimation of an examinee ability
 * under the two-parameter logistic item characteristic curve
 * model using LSAT-6 data.
 */
		public void modelBayesianItemParameters() {
			double[] jirp = new double[NITEM];
			jirp[0] = 1.0; // read canned data
			jirp[1] = 0.0;
			jirp[2] = 1.0;
			jirp[3] = 0.0;
			jirp[4] = 1.0;
			double[] b = { -3.283, -1.317, -0.285, -1.766, -2.858 };
			double[] a = { 0.849, 0.758, 0.870, 0.736, 0.730 };
			int trall = -1; // set traces
			int trbme = -1;
			String yn = JOptionPane.showInputDialog
					("Trace all? Y/N");
			if (yn.equalsIgnoreCase("y")) {
				trall = 0;
				trbme = 0;
			}
			else{
				String ynbme = JOptionPane.showInputDialog
						("Trace BME? Y/N");
				if (ynbme.equalsIgnoreCase("y")) {
					trbme = 0;
				}
			} // end if (yn ... else
			String smaxit = JOptionPane.showInputDialog
					("Enter maximum number of iterations to do. ");
			int maxit = Integer.parseInt(smaxit);
			System.out.println("All= " + trall +
					"  BME= " + trbme +
					"  Max Iterations= " + maxit);
			double theta = 0.0; // initial estimate of examinee's ability
			// estimate ability
			for (int nit = 1; nit <= maxit; nit++) {
				double sumnum = 0.0; // clear sums
				double sumdem = 0.0;
				System.out.println("Iteration= " + nit);
				for (int i = 0; i <= NITEM - 1; i++) { // item loop
					double uij = jirp[i];
					double dev = a[i] * (theta - b[i]);
					double phat = 1.0 / (1.0 + Math.exp(-dev));
					double wij = phat * (1.0 - phat);
					double vij = uij - phat;
					sumnum += a[i] * vij;
					sumdem += a[i] * a[i] * wij;
					if (trall == 0) {
						System.out.println("Item= " + (i + 1));
						System.out.println("Difficulty= " + b[i] +
								"  Discrimination= " + a[i] +
								"  Theta= " + theta);
						System.out.println("Uij= " + uij +
								"  Dev= " + dev +
								"  P Hat= " + phat);
						System.out.println("Wij= " + wij +
								"  Vij= " + vij);
						System.out.println("Numerator Sum= " + sumnum +
								"  Denominator Sum= " + sumdem);
					}
				} // end for (int i ... (BASIC 2190)
				sumnum -= (theta - HPM) / HPVAR;
				sumdem = - sumdem - (1.0 / HPVAR); // match equation 7.39
				if (trall == 0) {
					System.out.println("After adding proir terms");
					System.out.println("Numerator Sum= " + sumnum +
							"  Denominator Sum= " + sumdem);
				}
				double delta = sumnum / sumdem;
				if (trall == 0 || trbme == 0) {
					System.out.println("Change in Theta= " + delta);
				}
				if (Math.abs(delta) < BIGT) {
					// normal processing
				}
				else { // protect against big change in theta
					if (delta > 0.0) {
						delta = BIGT;
					}
					else {
						delta = - BIGT;
					} // end if (delta ...
				} // end if (Math.abs(delta) ...
				theta -= delta; // change subtracted as per equation 7.40
				if (trall == 0 || trbme == 0) {
					System.out.println("Theta= " + theta +
							"  Change= " + delta);
				} // end if (trall ...
				if (Math.abs(delta) < 0.05) {
					break; // exit for (int nit ...
				}
				if (nit == maxit) {
					System.out.println("Reached Max Iterations");
				}
			} // end for (int nit ... (BASIC 2460)
			System.out.println("Estimated Ability= " + theta);
			System.exit(0); // last line
		} // end public static ...

/**
 * A program to implement Bayesian expected a posteriori ability
 * estimation using LSAT-6 data.
 */

		public void bayesianAPrioriPersonParameter() {

			double[] u = {1.0, 0.0, 1.0, 0.0, 1.0};
			double[] x = { -4.000, -3.111, -2.222, -1.333, -0.4444,
					0.4444, 1.333, 2.222, 3.111, 4.000 };
			double[] ak = { 0.000119, 0.002805, 0.03002, 0.1458, 0.3213,
					0.3213, 0.1458, 0.03002, 0.002805, 0.000119 };
			double[] b = { -3.283, -1.317, -0.285, -1.766, -2.858 };
			double[] a = { 0.849, 0.758, 0.870, 0.736, 0.730 };
			double[] lxk = new double[NQ];
			int trall = -1; // set traces
			int trls = -1;
			String yn = JOptionPane.showInputDialog("Trace all? Y/N");
			if (yn.equalsIgnoreCase("y")) {
				trall = 0;
				trls = 0;
			} else {
				String ynls = JOptionPane.showInputDialog("Trace L and Sums? Y/N");
				if (ynls.equalsIgnoreCase("y")) {
					trls = 0;
				}
			} // end if (yn ... else

			System.out.println("All= " + trall + "  L and Sums= " + trls);

			for (int k = 0; k <= NQ - 1; k++) {
				lxk[k] = 1.0;
			} // end for (int k ...

			// Bayesian expected a posteriori estimation of ability
			for (int k = 0; k <= NQ - 1; k++) {
				for (int i = 0; i <= NITEM - 1; i++) {
					double dev = - a[i] * (x[k] - b[i]); // evaluate 2pm icc
					double p = 1.0 / (1.0 + Math.exp(dev));
					double q = 1.0 - p;
					double pitem = p;
					if (u[i] == 0.0) {
						pitem = q;
					}
					lxk[k] *= pitem;
					if (trall == 0) {
						System.out.println("K= " + (k + 1) + "  I= " + (i + 1) +
								"  U= " + u[i] + "  P= " + pitem + "  L(X)= " + lxk[k]);
					} // end if (trall ...
				} // end for (int i ... (BASIC 370)
				if (trall == 0 || trls == 0) {
					System.out.println("K= " + (k + 1) + "  L(X)= " + lxk[k]);
				}
			} // end for (int k ... (BASIC 390)
			double snum = 0.0;
			double sden = 0.0;
			for (int k = 0; k <= NQ - 1; k++) {
				double la = lxk[k]*ak[k];
				double sla = x[k] * la;
				snum += sla;
				sden += la;
				if (trall == 0) {
					System.out.println("K= " + (k + 1) + "  LXK(K)= " + lxk[k] +
							"  AK[K]= " + ak[k] + "  X(K)= " + x[k]);
					System.out.println("LA= " + la + "  XLA= " + sla);
				} // end if (trall ...
				if (trall == 0 || trls == 0) {
					System.out.println("K= " + (k + 1) + "  Sum Numerator= " + snum +
							"  Sum Denominator = " + sden);
				} // end if (trall ...
			} // end for (int k ... (BASIC 510)
			double theta = snum / sden;
			System.out.println("Estimated Theta= " + theta);
			System.exit(0); // last line
		} // end public static ...
	} // end public class ...
