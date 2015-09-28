package org.paces.Stata.IRTBaker;

import javax.swing.*;

/**
 * Created by billy on 9/27/15.
 */
public class MMLERasch {


/**
 * AppendixH.java
 * Estimation of Marginal Maximumm Likelihood Rasch
 * Inverse Whole matrix
 */
		static final int NI = 5; // Number of items I
		static final int NG = 6; // Number of score Groups
		static final int NL = 10; // Number of theta Levels
		static final int NIT = 24; // Number if ITerations 4*(NI+1)
		//static final int MNC = 25; // Maximum Number of em Cycles
		static final double CRITIT = 0.0001; // CRITerion of IT
		static final double CRITNC = 0.001; // CRITerion of NC
		//do not employ CRITNC
//g=0(1)5
		public MMLERasch(String[] args) {
			this.raschMMLE();
			this.raschAlternate();
		}

		public void raschMMLE() {
			int trall = -1; // set up trace of computations
			String yn = JOptionPane.showInputDialog
					("Trace all? Y/N");
			if (yn.equalsIgnoreCase("y")) {
				trall = 0;
			}
			else{
				// normal processing
			} // end if else
			String smaxit = JOptionPane.showInputDialog
					("Enter number of EM cycles to do. ");
			int MNC = Integer.parseInt(smaxit);
			double[] phix = new double[NL];
			double[][] cx = new double[NG][NL];
			double[] prodp = new double[NL];
			double[] scxphix = new double[NL];
			double[][] lx = new double[NG][NL];
			double[][] r1 = new double[NL][NI];
			double[][] r0 = new double[NL][NI];
			double[][] p = new double[NL][NI];
			double[] db = new double[NI];
			double[] ddb = new double[NI];
			double[] ddc = new double[NI]; // ddc is dadb
			double[][] ddi = new double[NI+1][NI+1]; // ddi is dd-inverse
			double[] deltab = new double[NI];
			double[] oldb = new double[NI];
			double[] x = {
					-4.5, -3.5, -2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 3.5, 4.5
			};
			double[][] rp = {
					{ 0.0, 0.0, 0.0, 0.0, 0.0 },
					{ 10.0, 1.0, 1.0, 2.0, 6.0 },
					{ 62.0, 24.0, 7.0, 28.0, 49.0 },
					{ 212.0, 109.0, 63.0, 139.0, 188.0 },
					{ 342.0, 277.0, 184.0, 296.0, 329.0 },
					{ 298.0, 298.0, 298.0, 298.0, 298.0 }
			};
			double[] n = {
					3.0, 20.0, 85.0, 237.0, 357.0, 298.0
			};
			// initial values
			double a = 1.0;
			double[] b = {
					0.0, 0.0, 0.0, 0.0, 0.0
			};
			double sumn = 0.0;
			for (int g = 0; g <= NG - 1; g++) {
				sumn += n[g];
			} // end for (int g ...
			for (int i = 0; i <= NI - 1; i++) {
				for (int g = 0; g <= NG - 1; g++) {
					b[i] += rp[g][i];
				} // end for (int g ...
				b[i] = Math.log((sumn - b[i])/ b[i]);
				System.out.println("I= " + (i + 1) +
						"  B[I]=  " + b[i]);
			} // end for (int i ...
			for (int l = 0; l <= NL - 1; l++) {
				phix[l] = (1.0 / Math.sqrt(2.0 * Math.PI)) *
						Math.exp(-(1.0 / 2.0) * Math.pow(x[l], 2));
			} // end for (int l ...
			for (int nc = 0; nc <= MNC - 1; nc++) {
				System.out.println("Cycle= " + (nc + 1));
				for (int i = 0; i <= NI - 1; i++) {
					oldb[i] = b[i];
				} // end for (int i ...
				double olda = a;
//    for CRITNC
				for (int l = 0; l <= NL - 1; l++) {
					prodp[l] = 1.0;
					for (int i = 0; i <= NI - 1; i++) {
						p[l][i] = 1.0 / (1.0 + Math.exp(-a * (x[l] - b[i])));
						prodp[l] *= p[l][i];
					} // end for (int i ...
				} // end for (int l ...
				for (int g = 0; g <= NG - 1; g++) {
					for (int l = 0; l <= NL - 1; l++) {
						cx[g][l] = prodp[l] * Math.exp(-(NI - g) * a * x[l]);
					} // end for (int l ...
				} // end for (int g ...
				for (int g = 0; g <= NG - 1; g++) {
					scxphix[g] = 0.0;
					for (int l = 0; l <= NL - 1; l++) {
						scxphix[g] += cx[g][l] * phix[l];
					} // end for (int l ...
				} // end for (int g ...
				for (int g = 0; g <= NG - 1; g++) {
					for (int l = 0; l <= NL - 1; l++) {
						lx[g][l] = cx[g][l] * phix[l] / scxphix[g];
					} // end for (int l ...
				} // end for (int g ...
				for (int l = 0; l <= NL - 1; l++) {
					for (int i = 0; i <= NI - 1; i++) {
						r1[l][i] = 0.0;
						r0[l][i] = 0.0;
						for (int g = 0; g <= NG - 1; g++) {
							r1[l][i] += rp[g][i] * lx[g][l];
							r0[l][i] += (n[g] - rp[g][i]) * lx[g][l];
						} // end for (int g ...
					} // end for (int l ...
				} // end for (int g ...
				for (int it = 0; it <= NIT - 1; it++) {
					for (int i = 0; i <= NI - 1; i++) {
						db[i] = 0.0;
						ddb[i] = 0.0;
						for (int l = 0; l <= NL - 1; l++) {
							p[l][i] = 1.0 / (1.0 + Math.exp(-a * (x[l] - b[i])));
							db[i] += a * (r0[l][i] * p[l][i] - r1[l][i] *
									(1.0 - p[l][i]));
							ddb[i] -= a * a * p[l][i] * (1.0 - p[l][i]) *
									(r0[l][i] + r1[l][i]);
						} // end for (int l ...
						if (trall == 0) {
							System.out.println("IT= " + (it + 1) +
									"  I= " + (i + 1) +
									"  B[I]= " + b[i] +
									"  DB[I]= " + db[i] +
									"  DDB[I]= " + ddb[i]);
						}
						else {
							// normal processing
						} // end if (trall ... else ...
					} // end for (int i ...
					double da = 0.0;
					double dda = 0.0;
					for (int l = 0; l <= NL - 1; l++) {
						for (int i = 0; i <= NI - 1; i++) {
							p[l][i] = 1.0 / (1.0 + Math.exp(-a * (x[l] - b[i])));
							da += (x[l] - b[i]) * (r1[l][i] * (1.0 - p[l][i]) -
									r0[l][i] * p[l][i]);
							dda -= (x[l] - b[i]) * (x[l] - b[i]) * p[l][i] *
									(1.0 - p[l][i]) * (r0[l][i] + r1[l][i]);
						} // end for (int i ...
					} // end for (int l ...
					if (trall == 0) {
						System.out.println("IT= " + (it + 1) +
								"  A= " + a +
								"  DA= " + da +
								"  DDA= " + dda);
					}
					else {
						// normal processing
					} // end if (trall ... else ...
					for (int i = 0; i <= NI - 1; i++) {
						ddc[i] = 0.0;
						for (int l = 0; l <= NL - 1; l++) {
							p[l][i] = 1.0 / (1.0 + Math.exp(-a * (x[l] - b[i])));
							ddc[i] += (r0[l][i] * p[l][i] - r1[l][i] *
									(1.0 - p[l][i])) + a * p[l][i] * (1.0 - p[l][i]) *
									(x[l] - b[i]) * (r0[l][i] + r1[l][i]);
						} // end for (int l ...
						if (trall == 0) {
							System.out.println("I= " + (i + 1) +
									"  DDC[I]= " + ddc[i]);
						}
						else {
							// normal processing
						} // end if (trall ... else ...
					} // end for (int i ..

					if (trall == 0) {
						System.out.println("Inverse Hessin");
					}
					else {
						// normal processing
					} // end if (trall ... else ...
					ddi[NI][NI] = 0.0;
					for (int i = 0; i <= NI - 1; i++) {
						ddi[NI][NI] += ddc[i] * (1.0 / ddb[i]) * ddc[i];
					} // end for (int i ...
					ddi[NI][NI] = 1.0 / (dda - ddi[NI][NI]);
					for (int i = 0; i <= NI - 1; i++) {
						ddi[i][NI] = -(1.0 / ddb[i]) * ddc[i] * ddi[NI][NI];
						ddi[NI][i] = -ddi[NI][NI] * ddc[i] / ddb[i];
					} // end for (int i ...
					for (int i = 0; i <= NI - 1; i++) {
						for (int ii = 0; ii <= NI - 1; ii++) {
							ddi[i][ii] = (1.0 / ddb[i]) * ddc[i] * ddi[NI][NI] *
									(ddc[ii] / ddb[ii]);
						} // end for (int ii ...
					} // end for (int i ...
					for (int i = 0; i <= NI - 1; i++) {
						ddi[i][i] = (1.0 / ddb[i]) + ddi[i][i];
					} // end for (int i ...
					if (trall == 0) {
						for (int i = 0; i <= NI; i++){
							for (int ii = 0; ii <= NI; ii++){
								System.out.println("I= " + (i + 1) +
										"  II= " + (ii + 1) +
										"  DDI[I][II]= " + ddi[i][ii]);
							} // end for (int ii ...
						} // end for (int i ...
					}
					else {
						// normal processing
					} // end if (trall ... else ...
					for (int i = 0; i <= NI - 1; i++) {
						deltab[i] = (ddi[i][0] * db[0] + ddi[i][1] * db[1] +
								ddi[i][2] * db[2] + ddi[i][3] * db[3] +
								ddi[i][4] * db[4] + ddi[i][5] * da);
						b[i] -= deltab[i];
					} // end for (int i ...
					double deltaa = (ddi[NI][0] * db[0] + ddi[NI][1] * db[1] +
							ddi[NI][2] * db[2] + ddi[NI][3] * db[3] +
							ddi[NI][4] * db[4] + ddi[NI][5] * da);
					a -= deltaa;
					if (trall == 0 ) {
						for (int i = 0; i <= NI - 1; i++) {
							System.out.println("I= " + (i + 1) +
									"  B[I]= " + b[i]);
						} // end for (int i ...
						System.out.println("A= " + a);
					}
					else {
						// normal processing
					} // end if (trall ... else ...
					double diffit = 0.0;
					for (int i = 0; i <= NI - 1; i++) {
						if (Math.abs(deltab[i]) > diffit) {
							diffit = Math.abs(deltab[i]);
						}
						else {
							// normal processing
						} // end if (Maht.abs ... else ...
					} // end for (int i ...
					if (Math.abs(deltaa) > diffit) {
						diffit = Math.abs(deltaa);
					}
					else {
						// normal processing
					} // end if (Math.abs ... else ...
					if (trall == 0 ) {
						System.out.println("IT= " + (it + 1) +
								"  DIFFIT= " + diffit);
					}
					else {
						// normal processing
					} // end if (trall ... else ...
					if (diffit < CRITIT) {
						break; // exit for (int it ...
					}
					else {
						// normal processing
					} // end if (diffit ... else ...
				} // end (int it ...
				if (trall == -1 ) {
					for (int i = 0; i <= NI - 1; i++) {
						System.out.println("I= " + (i + 1) +
								"  B[I]= " + b[i]);
					} // end for (int i ...
					System.out.println("A= " + a);
				}
				else {
					// normal processing
				} // end if (trall ... else ...
				double diffnc = Math.abs(a - olda);
				for (int i = 0; i <= NI - 1; i++) {
					if (Math.abs(b[i] - oldb[i]) > diffnc) {
						diffnc = Math.abs(b[i] - oldb[i]);
					}
					else {
						// normal processing
					} // end if (Math.abs(b[i] ... else ...
				} // end for (int i ...
				System.out.println("NC= " + (nc + 1) +
						"  DIFFNC=" + diffnc);
			} // end for (int nc ...
		} // end public static ...

/**
 * AppendixH1.java
 * Estimation of Marginal Maximumm Likelihood Rasch
 * Simpler Relaxation Solution
 */
			public void raschAlternate() {
				int trall = -1; // set up trace of computations
				String yn = JOptionPane.showInputDialog
						("Trace all? Y/N");
				if (yn.equalsIgnoreCase("y")) {
					trall = 0;
				}
				else{
					// normal processing
				} // end if else
				String smaxit = JOptionPane.showInputDialog
						("Enter number of EM cycles to do. ");
				int MNC = Integer.parseInt(smaxit);
				double[] phix = new double[NL];
				double[][] cx = new double[NG][NL];
				double[] prodp = new double[NL];
				double[] scxphix = new double[NL];
				double[][] lx = new double[NG][NL];
				double[][] r1 = new double[NL][NI];
				double[][] r0 = new double[NL][NI];
				double[][] p = new double[NL][NI];
				double[] db = new double[NI];
				double[] ddb = new double[NI];
				double[] oldb = new double[NI];
				double[] x = {
						-4.5, -3.5, -2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 3.5, 4.5
				};
				double[][] rp = {
						{ 0.0, 0.0, 0.0, 0.0, 0.0 },
						{ 10.0, 1.0, 1.0, 2.0, 6.0 },
						{ 62.0, 24.0, 7.0, 28.0, 49.0 },
						{ 212.0, 109.0, 63.0, 139.0, 188.0 },
						{ 342.0, 277.0, 184.0, 296.0, 329.0 },
						{ 298.0, 298.0, 298.0, 298.0, 298.0 }
				};
				double[] n = {
						3.0, 20.0, 85.0, 237.0, 357.0, 298.0
				};
				// initial values
				double a = 1.0;
				double[] b = {
						0.0, 0.0, 0.0, 0.0, 0.0
				};
				double sumn = 0.0;
				for (int g = 0; g <= NG - 1; g++) {
					sumn += n[g];
				} // end for (int g ...
				for (int i = 0; i <= NI - 1; i++) {
					for (int g = 0; g <= NG - 1; g++) {
						b[i] += rp[g][i];
					} // end for (int g ...
					b[i] = Math.log((sumn - b[i])/ b[i]);
					System.out.println("I= " + (i + 1) +
							"  B[I]=  " + b[i]);
				} // end for (int i ...
				for (int l = 0; l <= NL - 1; l++) {
					phix[l] = (1.0 / Math.sqrt(2.0 * Math.PI)) *
							Math.exp(-(1.0 / 2.0) * Math.pow(x[l], 2));
				} // end for (int l ...
				for (int nc = 0; nc <= MNC - 1; nc++) {
					System.out.println("Cycle= " + (nc + 1));
					for (int i = 0; i <= NI - 1; i++) {
						oldb[i] = b[i];
					} // end for (int i ...
					double olda = a;
//    for CRITNC
					for (int l = 0; l <= NL - 1; l++) {
						prodp[l] = 1.0;
						for (int i = 0; i <= NI - 1; i++) {
							p[l][i] = 1.0 / (1.0 + Math.exp(-a * (x[l] - b[i])));
							prodp[l] *= p[l][i];
						} // end for (int i ...
					} // end for (int l ...
					for (int g = 0; g <= NG - 1; g++) {
						for (int l = 0; l <= NL - 1; l++) {
							cx[g][l] = prodp[l] * Math.exp(-(NI - g) * a * x[l]);
						} // end for (int l ...
					} // end for (int g ...
					for (int g = 0; g <= NG - 1; g++) {
						scxphix[g] = 0.0;
						for (int l = 0; l <= NL - 1; l++) {
							scxphix[g] += cx[g][l] * phix[l];
						} // end for (int l ...
					} // end for (int g ...
					for (int g = 0; g <= NG - 1; g++) {
						for (int l = 0; l <= NL - 1; l++) {
							lx[g][l] = cx[g][l] * phix[l] / scxphix[g];
						} // end for (int l ...
					} // end for (int g ...
					for (int l = 0; l <= NL - 1; l++) {
						for (int i = 0; i <= NI - 1; i++) {
							r1[l][i] = 0.0;
							r0[l][i] = 0.0;
							for (int g = 0; g <= NG - 1; g++) {
								r1[l][i] += rp[g][i] * lx[g][l];
								r0[l][i] += (n[g] - rp[g][i]) * lx[g][l];
							} // end for (int g ...
						} // end for (int l ...
					} // end for (int g ...
					for (int i = 0; i <= NI - 1; i++) {
						db[i] = 0.0;
						ddb[i] = 0.0;
						for (int it = 0; it <= NIT - 1; it++) {
							for (int l = 0; l <= NL - 1; l++) {
								p[l][i] = 1.0 / (1.0 + Math.exp(-a * (x[l] - b[i])));
								db[i] += a * (r0[l][i] * p[l][i] - r1[l][i] *
										(1.0 - p[l][i]));
								ddb[i] -= a * a * p[l][i] * (1.0 - p[l][i]) *
										(r0[l][i] + r1[l][i]);
							} // end for (int l ...
							if (trall == 0) {
								System.out.println("IT= " + (it + 1) +
										"  I= " + (i + 1) +
										"  B[I]= " + b[i] +
										"  DB[I]= " + db[i] +
										"  DDB[I]= " + ddb[i]);
							}
							else {
								// normal processing
							} // end if (trall ... else ...
							b[i] -= db[i] / ddb[i];
							if (Math.abs(db[i] / ddb[i]) < CRITIT) {
								break; // exit for (int it ...
							}
							else {
								// normal processing
							} // end if (Math.ab ...
						} // end for (int it ...
					} // end for (int i ...
					for (int it = 0; it <= NIT - 1; it++) {
						double da = 0.0;
						double dda = 0.0;
						for (int l = 0; l <= NL - 1; l++) {
							for (int i = 0; i <= NI - 1; i++) {
								p[l][i] = 1.0 / (1.0 + Math.exp(-a * (x[l] - b[i])));
								da += (x[l] - b[i]) * (r1[l][i] * (1.0 - p[l][i]) -
										r0[l][i] * p[l][i]);
								dda -= (x[l] - b[i]) * (x[l] - b[i]) * p[l][i] *
										(1.0 - p[l][i]) * (r0[l][i] + r1[l][i]);
							} // end for (int i ...
						} // end for (int l ...
						if (trall == 0) {
							System.out.println("IT= " + (it + 1) +
									"  A= " + a +
									"  DA= " + da +
									"  DDA= " + dda);
						}
						else {
							// normal processing
						} // end if (trall ... else ...
						a -= da / dda;
						if (Math.abs(da / dda) < CRITIT) {
							break; // exit for (int it ...
						}
						else {
							// normal processing
						} // end if (Math.ab ...
					} // end for (int it ...
					for (int i = 0; i <= NI - 1; i++) {
						System.out.println("I= " + (i + 1) +
								"  B[I]= " + b[i]);
					} // end for (int i ...
					System.out.println("A= " + a);
					double diffnc = Math.abs(a - olda);
					for (int i = 0; i <= NI - 1; i++) {
						if (Math.abs(b[i] - oldb[i]) > diffnc) {
							diffnc = Math.abs(b[i] - oldb[i]);
						}
						else {
							// normal processing
						} // end if (Math.abs(b[i] ... else ...
					} // end for (int i ...
					System.out.println("NC= " + (nc + 1) +
							"  DIFFNC=" + diffnc);
				} // end for (int nc ...
			} // end public static ...
		} // end public class ...
