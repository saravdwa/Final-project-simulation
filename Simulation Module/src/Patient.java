import static java.lang.StrictMath.max;

public class Patient {
    // object class
    public int nr;
    public int patientType;    // (0=none, 1=elective, 2=urgent)
    public int scanType;       // elective: (0=none), urgent: (0=brain, 1=lumbar, 2=cervival, 3=abdomen, 4=others)
    public int callWeek;       // week of arrival (elective: call, urgent: actual arrival)
    public int callDay;        // day of arrival (elective: call, urgent: actual arrival)
    public double callTime;    // time of arrival (elective: call, urgent: actual arrival) (in hours)
    public int scanWeek;       // week of appointment
    public int scanDay;        // day of appointment
    public int slotNr;         // slot number of appointment
    public double appTime;     // time of appointment (elective: according to rule, urgent: slot start time) (in hours)
    public double tardiness;   // (in hours)
    public boolean isNoShow;
    public double scanTime;    // actual start time of the scan
    public double duration;    // actual duration of the scan

    public Patient(int nr_, int patientType_, int scanType_, int callWeek_, int callDay_, double callTime_, double tardiness_, boolean isNoShow_, double duration_){
        nr = nr_;
        patientType = patientType_;
        scanType = scanType_;
        callWeek = callWeek_;
        callDay = callDay_;
        callTime = callTime_;
        tardiness = tardiness_;
        isNoShow = isNoShow_;
        duration = duration_;

        //unplanned
        scanWeek = -1;
        scanDay = -1;
        slotNr = -1;
        appTime = -1;
        scanTime = -1.0;
    }

    public double getAppWT(){
        if(slotNr != -1){
            return (double)(((scanWeek-callWeek)*7 + scanDay - callDay)*24 + appTime - callTime); // in hours
        }else{
            System.out.println("CAN NOT CALCULATE APPOINTMENT WT OF PATIENT %d" + nr);
            System.exit(1);
        }
        return 0;
    }

    public double getScanWT(){
        if(scanTime != 0){
            double wt = 0;
            if(patientType == 1){ // elective
                wt = scanTime - (appTime + tardiness);
            }else{ // urgent
                wt = scanTime - callTime;
            }
            return max(0.0,wt);
        }else{
            System.out.println("CAN NOT CALCULATE SCAN WT OF PATIENT %d" + nr);  // in hours
            System.exit(1);
        }
        return 0;
    }

}

