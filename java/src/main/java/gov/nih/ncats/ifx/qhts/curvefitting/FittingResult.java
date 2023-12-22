package gov.nih.ncats.ifx.qhts.curvefitting;
 
import java.io.Serializable;
  
public class FittingResult implements Serializable{

	private static final long serialVersionUID = 1L;
	
	
	private Integer maskNo;	
	private String maskFlag;
	
	private Double zeroActivity;
	private Double infActivity;
	private Double logAc50;
	private Double hillCoef;	
	
	private Double curveClass2;
	private Double ac50;
	private Double efficacy;
	private Double maxResponse;
	private Double pValue;
	private Double r2;
	
	
	public FittingResult() {
		super();
	}

	/*
	@Override
	public String toString() {
		return "FittingResult [\nac50=" + ac50 + ",\n curveClass2=" + curveClass2 + ",\n efficacy=" + efficacy + ",\n hillCoef="
				+ hillCoef + ",\n infActivity=" + infActivity + ",\n logAc50=" + logAc50 + ",\n maskFlag=" + maskFlag
				+ ",\n maskNo=" + maskNo + ",\n maxResponse=" + maxResponse + ",\n pValue=" + pValue + ",\n r2=" + r2
				+ ",\n zeroActivity=" + zeroActivity + "\n]";
	}
	*/
	
	@Override
	public String toString() {
		return ac50 + "\n" + curveClass2 + "\n" + efficacy + "\n" +
			   hillCoef + "\n" + infActivity + "\n" + logAc50 + "\n" + 
			   maskFlag + "\n" + maskNo + "\n" + maxResponse + "\n" + 
			   pValue + "\n" + r2 + "\n" + zeroActivity + "\n";
	}


	public Integer getMaskNo() {
		return maskNo;
	}


	public void setMaskNo(Integer maskNo) {
		this.maskNo = maskNo;
	}


	public String getMaskFlag() {
		return maskFlag;
	}


	public void setMaskFlag(String maskFlag) {
		this.maskFlag = maskFlag;
	}


	public Double getZeroActivity() {
		return zeroActivity;
	}


	public void setZeroActivity(Double zeroActivity) {
		this.zeroActivity = zeroActivity;
	}


	public Double getInfActivity() {
		return infActivity;
	}


	public void setInfActivity(Double infActivity) {
		this.infActivity = infActivity;
	}


	public Double getLogAc50() {
		return logAc50;
	}


	public void setLogAc50(Double logAc50) {
		this.logAc50 = logAc50;
	}


	public Double getHillCoef() {
		return hillCoef;
	}


	public void setHillCoef(Double hillCoef) {
		this.hillCoef = hillCoef;
	}


	public Double getCurveClass2() {
		return curveClass2;
	}


	public void setCurveClass2(Double curveClass2) {
		this.curveClass2 = curveClass2;
	}


	public Double getAc50() {
		return ac50;
	}


	public void setAc50(Double ac50) {
		this.ac50 = ac50;
	}


	public Double getEfficacy() {
		return efficacy;
	}


	public void setEfficacy(Double efficacy) {
		this.efficacy = efficacy;
	}


	public Double getMaxResponse() {
		return maxResponse;
	}


	public void setMaxResponse(Double maxResponse) {
		this.maxResponse = maxResponse;
	}


	public Double getpValue() {
		return pValue;
	}


	public void setpValue(Double pValue) {
		this.pValue = pValue;
	}


	public Double getR2() {
		return r2;
	}


	public void setR2(Double r2) {
		this.r2 = r2;
	}

}
