package org.paces.Stata.IRTBaker;

import java.util.HashMap;
import java.util.Map;

/**
 * Created by billy on 9/28/15.
 */
public class BakerDemos {

	public static void main(String[] args) {

		Map<String, Integer> argMap = new HashMap<String, Integer>();
		argMap.put("bayes", 1);
		argMap.put("gibbs", 2);
		argMap.put("grm", 3);
		argMap.put("jmle", 4);
		argMap.put("mixed", 5);
		argMap.put("mlei", 6);
		argMap.put("nominal", 7);
		argMap.put("mlep", 8);
		argMap.put("mmlei", 9);
		argMap.put("mmler", 10);
		argMap.put("multi", 11);
		argMap.put("", 4);


		switch(argMap.get(args[0])) {
			case 1 :
				BakerDemos.bayes(args);
				break;
			case 2 :
				BakerDemos.gibbs(args);
				break;
			case 3 :
				BakerDemos.gradedResponseModel(args);
				break;
			case 4 :
				BakerDemos.jmleRasch(args);
				break;
			case 5 :
				BakerDemos.mixed(args);
				break;
			case 6 :
				BakerDemos.mleItems(args);
				break;
			case 7 :
				BakerDemos.mleNominal(args);
				break;
			case 8 :
				BakerDemos.mlePerson(args);
				break;
			case 9 :
				BakerDemos.mmleItems(args);
				break;
			case 10 :
				BakerDemos.mmleRasch(args);
				break;
			case 11 :
				BakerDemos.multiGroup(args);
				break;
			default :
				BakerDemos.jmleRasch(args);
				break;
		}

	}

	public static void bayes(String[] args) {
		BayesModelItemParameters bayesmod = new BayesModelItemParameters(args);
	}
	public static void gibbs(String[] args) {
		GibbsSampler gibbsMCMC = new GibbsSampler(args);
	}
	public static void gradedResponseModel(String[] args) {
		GRMItemParameters grmParams = new GRMItemParameters(args);
	}
	public static void jmleRasch(String[] args){
		JMLERasch jmle = new JMLERasch(args);
	}
	public static void mixed(String[] args){
		MixedModels jmle = new MixedModels(args);
	}
	public static void mleItems(String[] args){
		MLEItemParameters jmle = new MLEItemParameters(args);
	}
	public static void mleNominal(String[] args){
		MLENominalResponse jmle = new MLENominalResponse(args);
	}
	public static void mlePerson(String[] args){
		MLEPersonParameter mmle = new MLEPersonParameter(args);
	}
	public static void mmleItems(String[] args){
		MMLEItemParameters mmleItems = new MMLEItemParameters(args);
	}
	public static void mmleRasch(String[] args){
		MMLERasch mmle = new MMLERasch(args);
	}
	public static void multiGroup(String[] args){
		MultipleGroups multi = new MultipleGroups(args);
	}
}
