import java.io.FileInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import static java.lang.Double.max;

public class Simulation {

    // Variables and parameters
    String inputFileName;
    int D = 6;                         // number of days per week (NOTE: Sunday not included! so do NOT use to calculate appointment waiting time)
    int amountOTSlotsPerDay = 10;      // number of overtime slots per day
    int S = 32 + amountOTSlotsPerDay;  // number of slots per day
    double slotLenght = 15.0 / 60.0;             // duration of a slot (in hours)
    double lambdaElective = 28.345;
    double meanTardiness = 0;
    double stdevTardiness = 2.5;
    double probNoShow = 0.02;
    double meanElectiveDuration = 15;
    double stdevElectiveDuration = 3;
    double[] lambdaUrgent = new double[]{2.5, 1.25};
    double[] probUrgentType = new double[]{0.7, 0.1, 0.1, 0.05, 0.05};
    double[] cumulativeProbUrgentType = new double[]{0.7, 0.8, 0.9, 0.95, 1.0};
    double[] meanUrgentDuration  = new double[] {15, 17.5, 22.5, 30, 30};
    double[] stdevUrgentDuration  = new double[]{2.5, 1, 2.5, 1, 4.5};
    double weightEl = 1.0 / 168.0;      // objective weight elective appointment wait time
    double weightUr = 1.0 / 9.0;        // objective weight urgent scan wait time

    int W, R;                          // number of weeks 😊 runs lenght) and number of replications (set their values yourself in the initalization method!)
    int d,s,w,r;
    int rule;                          // the appointment scheduling rule
    Slot[][] weekSchedule = new Slot[D][s];  // array of the cyclic slot schedule (days-slots)
    // TODO: moet nog GEINITEERD WORDEN?


    // Variables within one simuation
    List<Patient> patients = new ArrayList<>(); // patient list
    double[] movingAvgElectiveAppWT;    // moving average elective appointment waiting time
    double[] movingAvgElectiveScanWT;   // moving average elective scan waiting time
    double[] movingAvgUrgentScanWT;     // moving average urgent scan waiting time
    double[] movingAvgOT;               // moving average overtime
    double avgElectiveAppWT;           // average elective appointment waiting time
    double avgElectiveScanWT;          // average elective scan waiting time
    double avgUrgentScanWT;            // average urgent scan waiting time
    double avgOT;                      // average overtime
    int numberOfElectivePatientsPlanned;
    int numberOfUrgentPatientsPlanned;

    // Initialization of a "simulation" object --> dit zijn de 2 constructors denk ik (een gevulde en een lege)
        public Simulation(int D, int S, int W){
        // Set test case variables
        //TODO: set these variables to the correct values
        inputFileName = "/Users/tinemeersman/Documents/project SMA 2022 student code /input-S1-14.txt";  // input file with schedule
        W = 10;                      // number of weeks to simulate = run lenght
        R = 1;                      // number of replications
        rule = 1;                   // the appointment scheduling rule to apply

        // Initialize variables
        avgElectiveAppWT = 0;
        avgElectiveScanWT = 0;
        avgUrgentScanWT = 0;
        avgOT = 0;
        numberOfElectivePatientsPlanned = 0;
        numberOfUrgentPatientsPlanned = 0;

        // Initialize arrays
            weekSchedule = new Slot[D];
        for(this.d = 0; this.d < D; this.d++){
            weekSchedule[this.d] = new Slot[S];
        }
        movingAvgElectiveAppWT = new double[W];
        movingAvgElectiveScanWT = new double[W];
        movingAvgUrgentScanWT = new double[W];
        movingAvgOT = new double[W];
    }

    public Simulation(){
    }


    public void runSimulations(){
        double electiveAppWT = 0;
        double electiveScanWT = 0;
        double urgentScanWT = 0;
        double OT = 0;
        double OV = 0;
        setWeekSchedule();
        // set cyclic slot schedule based on given input file
        Random r= new Random();
        System.out.printf("r \t elAppWT \t elScanWT \t urScanWT \t OT \t OV \n");
        // run R replications
        for(s = 0; s < R; s++){
            resetSystem();          // reset all variables related to 1 replication
            r.setSeed(s); // set seed value for random value generator
            runOneSimulation();     // run 1 simulation / replication
            electiveAppWT += avgElectiveAppWT;
            electiveScanWT += avgElectiveScanWT;
            urgentScanWT += avgUrgentScanWT;
            OT += avgOT;
            OV += avgElectiveAppWT / weightEl + avgUrgentScanWT / weightUr;
            System.out.printf("%d \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \n", r, avgElectiveAppWT, avgElectiveScanWT, avgUrgentScanWT, avgOT, avgElectiveAppWT / weightEl + avgUrgentScanWT / weightUr);
        }
        electiveAppWT = electiveAppWT / R;
        electiveScanWT = electiveScanWT / R;
        urgentScanWT = urgentScanWT / R;
        OT = OT / R;
        OV = OV / R;
        double objectiveValue = electiveAppWT / weightEl + urgentScanWT / weightUr;
        System.out.printf("Avg.: \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \n", electiveAppWT, electiveScanWT, urgentScanWT, OT, objectiveValue);

        // print results
        //FILE *file = fopen("/Users/tinemeersman/Documents/project SMA 2022 student code /output.txt", "a"); // TODO: use your own directory
        // TODO: print the output you need to a .txt file
        //fclose(file);
    }

    public void setWeekSchedule(){
        // Read and set the slot types (0=none, 1=elective, 2=urgent within normal working hours)

        InputStream input = null;
        try {
            input = new FileInputStream(inputFileName);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
        Scanner inputFile = new Scanner(input);
        int elementInt;
        while (inputFile.hasNext()){
            for(s = 0; s < 32; s++){
                for(d = 0; d < D; d++){
                    elementInt = inputFile.nextInt();
                    weekSchedule[d][s].slotType = elementInt; //moet nog geiniteerd worden dus werkt nog niet
                    weekSchedule[d][s].patientType = elementInt;
                }
            }

        }

        inputFile.close();

        // Set the type of the overtime slots (3=urgent in overtime)
        for(d = 0; d < D; d++){
            for(s = 32; s < S; s++){
                weekSchedule[d][s].slotType = 3;
                weekSchedule[d][s].patientType = 2;
            }
        }

        // set start and appoitnment time
        double time;
        for(d = 0; d < D; d++){
            time = 8; // start time slot schedule
            for(s = 0; s < S; s++){
                // start time slot
                weekSchedule[d][s].startTime = time;

                // appointment time slot
                if(weekSchedule[d][s].slotType != 1){    // all slot types not elective : appointment time = slot start time
                    weekSchedule[d][s].appTime = time;
                }else{                                   // elective slots: appointment time is set according to rule !
                    if(rule == 1){ // FIFO
                        weekSchedule[d][s].appTime = time;
                    }else if(rule == 2){
                        // TODO: Bailey-Welch rule
                    }else if(rule == 3){
                        // TODO: Blocking rule
                    }else if(rule == 4){
                        // TODO: Benchmark rule
                    }
                }

                //update time variable
                time += slotLenght;
                if(time == 12){ time = 13;} // skip to the end of the luchbreak
            }
        }
    }

    public void resetSystem(){
        // reset all variables related to 1 replication
        patients.clear();
        avgElectiveAppWT = 0;
        avgElectiveScanWT = 0;
        avgUrgentScanWT = 0;
        avgOT = 0;
        numberOfElectivePatientsPlanned = 0;
        numberOfUrgentPatientsPlanned = 0;

        for(w = 0; w < W; w++){
            movingAvgElectiveAppWT[w] = 0;
            movingAvgElectiveScanWT[w] = 0;
            movingAvgUrgentScanWT[w] = 0;
            movingAvgOT[w] = 0;
        }
    }

    void runOneSimulation(){
        generatePatients();     // create patient arrival events (elective patients call, urgent patient arrive at the hospital)
        schedulePatients();     // schedule urgent and elective patients in slots based on their arrival events => detrmine the appointment wait time
        sortPatientsOnAppTime();   // sort patients on their appointment time (unscheduled patients are grouped at the end of the list)

        // determine scan wait time per patient and overtime per day
        int prevWeek = 0; int prevDay = -1;
        int[] numberOfPatientsWeek = new int[] {0,0};
        int[] numberOfPatients= new int[] {0,0};
        double arrivalTime, wt;
        double prevScanEndTime = 0;
        boolean prevIsNoShow = false;
        // go over arrival events (i.e. the moment the patient arrives at the hospital)
        for (int i = 0; i < patients.size(); i++) {
            {
                Patient nextPatient = new patients[i];
                if(patients.scanWeek == -1){ // stop at the first unplanned patient
                    break;
                }
        }

            arrivalTime = (double) patient->appTime + patient->tardiness;
            // SCAN WT
            if(!patient->isNoShow){
                if(patient->scanWeek != prevWeek || patient->scanDay != prevDay){
                    patient->scanTime = arrivalTime;
                } else{
                    if(prevIsNoShow){
                        patient->scanTime =  max(weekSchedule[patient->scanDay][patient->slotNr].startTime, max(prevScanEndTime,arrivalTime)); // note we assume we wait at least 15minutes on a no-show patient to see whether he shows or is just late
                    }else{
                        patient->scanTime = max(prevScanEndTime,arrivalTime);
                    }
                }
                wt = patient->getScanWT();
                if(patient->patientType == 1){
                    movingAvgElectiveScanWT[patient->scanWeek] += wt;
                }else{
                    movingAvgUrgentScanWT[patient->scanWeek] += wt;
                }
                numberOfPatientsWeek[patient->patientType - 1]++;
                if(patient->patientType == 1){
                    avgElectiveScanWT += wt;
                }else{
                    avgUrgentScanWT += wt;
                }
                numberOfPatients[patient->patientType - 1]++;
            }

            // OVERTIME
            if(prevDay > -1 && prevDay != patient->scanDay){
                if(d == 3 || d == 5){
                    movingAvgOT[prevWeek] += max(0.0, prevScanEndTime - 13);
                }else{
                    movingAvgOT[prevWeek] += max(0.0, prevScanEndTime - 17);
                }
                if(d == 3 || d == 5){
                    avgOT += max(0.0, prevScanEndTime - 13);
                }else{
                    avgOT += max(0.0, prevScanEndTime - 17);
                }
            }

            // update moving averages if week ends
            if(prevWeek != patient->scanWeek){
                movingAvgElectiveScanWT[prevWeek] = movingAvgElectiveScanWT[prevWeek] / numberOfPatientsWeek[0];
                movingAvgUrgentScanWT[prevWeek] = movingAvgUrgentScanWT[prevWeek] / numberOfPatientsWeek[1];
                movingAvgOT[prevWeek] = movingAvgOT[prevWeek] / D;
                numberOfPatientsWeek[0] = 0;
                numberOfPatientsWeek[1] = 0;
            }

            //set prev patient
            if(patient->isNoShow){
                //prevScanEndTime stays the same, it is the end time of the patient before the no-show patient
                prevIsNoShow = true;
            }else{
                prevScanEndTime = patient->scanTime + patient->duration;
                prevIsNoShow = false;
            }
            prevWeek = patient->scanWeek;
            prevDay = patient->scanDay;
        }
        // update moving averages of the last week
        movingAvgElectiveScanWT[W-1] = movingAvgElectiveScanWT[W-1] / numberOfPatientsWeek[0];
        movingAvgUrgentScanWT[W-1] = movingAvgUrgentScanWT[W-1] / numberOfPatientsWeek[1];
        movingAvgOT[W-1] = movingAvgOT[W-1] / D;

        // calculate objective values
        avgElectiveScanWT = avgElectiveScanWT / numberOfPatients[0];
        avgUrgentScanWT = avgUrgentScanWT / numberOfPatients[1];
        avgOT = avgOT / (D * W);


        // print moving avg
    /*FILE *file = fopen("/Users/tinemeersman/Documents/project SMA 2022 student code /output-movingAvg.txt", "a"); // TODO: use your own directory
    fprintf(file,"week \t elAppWT \t elScanWT \t urScanWT \t OT \n");
    for(w = 0; w < W; w++){
        fprintf(file, "%d \t %.2f \t %.2f \t %.2f \t %.2f \n", w, movingAvgElectiveAppWT[w], movingAvgElectiveScanWT[w], movingAvgUrgentScanWT[w], movingAvgOT[w]);
    }
    fclose(file);*/
    }


    public int getRandomScanType(){

        Random s= new Random();

        float r;
        r = (float) s.nextInt(1000-1)/ 1000;

        int type = -1;
        for(int i = 0; i < 5 && type == -1; i++){
            if(r < cumulativeProbUrgentType[i]){ type = i; }
        }
        return type;
    }

    public void generatePatients(){
        double arrivalTimeNext;
        int counter = 0; // total number of patients so far
        int patientType, scanType, endTime;
        double callTime, tardiness, duration, lambda;
        boolean noShow;
        Random r = null;
        for(w=0; w < W; w++){
            for(d = 0; d < D; d++){ // not on Sunday
                // generate ELECTIVE patients for this day
                if(d < D-1){  // not on Saturday either
                    arrivalTimeNext = 8 + Exponential_distribution(lambdaElective,r) * (17-8);
                    while(arrivalTimeNext < 17){ // desk open from 8h until 17h
                        patientType = 1;                // elective
                        scanType = 0;                   // no scan type
                        callTime = arrivalTimeNext;     // set call time, i.e. arrival event time
                        tardiness = Normal_distribution(meanTardiness, stdevTardiness,r) / 60.0;       // in practice this is not known yet at time of call
                        noShow = Bernouilli_distribution(probNoShow,r);                                // in practice this is not known yet at time of call
                        duration = Normal_distribution(meanElectiveDuration, stdevElectiveDuration,r) / 60.0; // in practice this is not known yet at time of call
                        Patient patient = new Patient(counter,patientType,scanType, w, d, callTime, tardiness, noShow, duration);
                        patients.add(patient);
                        counter++;
                        arrivalTimeNext = arrivalTimeNext + Exponential_distribution(lambdaElective,r) * (17-8); // arrival time of next patient (if < 17h)
                    }
                }

                // generate URGENT patients for this day
                if(d == 3 || d == 5){
                    lambda = lambdaUrgent[1]; // on Wed and Sat, only half a day!
                    endTime = 12;
                }else{
                    lambda = lambdaUrgent[0];
                    endTime = 17;
                }
                arrivalTimeNext = 8 + Exponential_distribution(lambda, r) * (endTime-8);
                while(arrivalTimeNext < endTime){ // desk open from 8h until 17h
                    patientType = 2;                // urgent
                    scanType = getRandomScanType(); // set scan type
                    callTime = arrivalTimeNext;     // set arrival time, i.e. arrival event time
                    tardiness = 0;                  // urgent patients have an arrival time = arrival event time
                    noShow = false;                 // urgent patients are never no-show
                    duration = Normal_distribution(meanUrgentDuration[scanType], stdevUrgentDuration[scanType],r) / 60.0; // in practice this is not known yet at time of arrival
                    Patient patient = new Patient(counter, patientType, scanType, w, d, callTime, tardiness, noShow, duration);
                    //ArrayList<patient> pa = new ArrayList<patient>();
                    patients.add(patient);
                    counter++;
                    arrivalTimeNext = arrivalTimeNext + Exponential_distribution(lambda, r) * (endTime-8); // arrival time of next patient (if < 17h)
                }
            }
        }
    }

    private boolean Bernouilli_distribution(double probNoShow, Random r) {
        double j1 = (float) r.nextInt(1000+1)/1000;
        boolean answer = false;
        if (j1 < probNoShow)
            return answer; // false because not a no show
        else {
            answer = true;
            return answer; // true because a no show
        }
    }

    private double Normal_distribution(double meanTardiness, double stdevTardiness, Random r) {
        double v1, v2, t;
        do{
            v1 = (float) r.nextInt(1000+1)*2;
            v1 /= 1000;
            v1 -= 1;
            v2 = (float) r.nextInt(1000+1)*2;
            v2 /= 1000;
            v2 -= 1;
            t=v1*v1+v2*v2;
        }
        while(t>=1||t==0);
        double multiplier = Math.sqrt(-2*Math.log(t)/t);
        return (int) (v1 * multiplier * stdevTardiness + meanTardiness);

    }

    private double Exponential_distribution(double lambdaElective, Random r) {
        double j1 = (float) r.nextInt(1000+1)/1000;
        if (j1 == 0)
            j1 += 0.0001;
        double j2 = -Math.log(j1)/lambdaElective;

        return j2;
    }

    public int getNextSlotNrFromTime(int day, int patientType, double time){
        bool found = false;
        int slotNr = -1;
        for(s = 0; !found && s < S; s++){
            if(weekSchedule[day][s].appTime > time && patientType == weekSchedule[day][s].patientType){
                found = true;
                slotNr = s;
            }
        }
        if(!found){
            printf("NO SLOT EXISTS DURING TIME %.2f \n", time);
            exit(0);
        }
        return slotNr;
    }

    public void schedulePatients(){
        //sort arrival events (= patient list) on arrival time (call time for elective patients, arrival time for urgent)
        patients.sort([](const Patient &patient1, const Patient &patient2){
            if (patient1.callWeek != patient2.callWeek)
                return patient1.callWeek < patient2.callWeek;
            if (patient1.callDay != patient2.callDay)
                return patient1.callDay < patient2.callDay;
            if (patient1.callTime != patient2.callTime)
                return patient1.callTime < patient2.callTime;
            if (patient1.scanType == 2)                             // if arrival time same, urgent patient before elective patient
                return true;
            if (patient2.scanType == 2)
                return false;
            return true;
        });

        int week[2] = {0,0}; // week of the next available slot {elective,urgent}
        int day[2] = {0,0}; // day of the next available slot {elective,urgent}
        int slot[2] = {0,0}; // slotNr of the next available slot {elective,urgent}

        //find first slot of each patient type (note, we assume each day (i.e. also day 0) has at least one slot of each patient type!)
        //elective
        d = 0;
        bool found = false;
        for(s = 0; s < S && !found; s++){
            if(weekSchedule[d][s].patientType == 1){
                day[0] = d;
                slot[0] = s;
                found = true;
            }
        }
        //urgent
        found = false;
        for(s = 0; s < S && !found; s++){
            if(weekSchedule[d][s].patientType == 2){
                day[1] = d;
                slot[1] = s;
                found = true;
            }
        }

        // go over SORTED patient list and assign slots
        int previousWeek = 0; int numberOfElective = 0; int numberOfElectivePerWeek = 0;   // keep track of week to know when to update moving average elective appointment waiting time
        double wt; int slotNr;
        for(patient = patients.begin(); patient != patients.end(); patient++){
            //Patient *pat = &*patient;

            //set index i dependant on patient type
            int i = patient->patientType - 1;

            // if still within the planning horizon:
            if(week[i] < W){

                // determine week where we start searching for a slot
                if(patient->callWeek > week[i]){
                    week[i] = patient->callWeek;
                    day[i] = 0;
                    slot[i] = getNextSlotNrFromTime(day[i], patient->patientType, 0);           // note we assume there is at least one slot of each patient type per day => this line will find first slot of this type
                }
                // determine day where we start searching for a slot
                if(patient->callWeek == week[i] && patient->callDay > day[i]){
                    day[i] = patient->callDay;
                    slot[i] = getNextSlotNrFromTime(day[i], patient->patientType, 0);           // note we assume there is at least one slot of each patient type per day => this line will find first slot of this type
                }
                // determine slot
                if(patient->callWeek == week[i] && patient->callDay == day[i] && patient->callTime >= weekSchedule[day[i]][slot[i]].appTime){
                    // find last slot on day "day[i]"
                    found = false; slotNr = -1;
                    for(s = S - 1; s >= 0 && !found; s--){
                        if(weekSchedule[day[i]][s].patientType == patient->patientType){
                            found = true;
                            slotNr = s;
                        }
                    }
                    // urgent patients have to be treated on the same day either in normal hours or in overtime (!! make sure there are enough overtime slots)
                    // for elective patients: check if the patient call time is before the last slot, i.e. if the patient can be planned on day "day[i]"
                    if(patient->patientType == 2 || patient->callTime < weekSchedule[day[i]][slotNr].appTime){
                        slot[i] = getNextSlotNrFromTime(day[i], patient->patientType, patient->callTime);   // find the first elective slot after the call time on day "day[i]"
                    }else{
                        // determine the next day
                        if(day[i] < D - 1){
                            day[i] = day[i] + 1;
                        }else{
                            day[i] = 0;
                            week[i] = week[i] + 1;
                        }
                        if(week[i] < W){   // find the first slot on the next day (if within the planning horizon)
                            slot[i] = getNextSlotNrFromTime(day[i], patient->patientType, 0);
                        }
                    }
                }

                // schedule the patient
                patient->scanWeek = week[i];
                patient->scanDay = day[i];
                patient->slotNr = slot[i];
                patient->appTime = weekSchedule[day[i]][slot[i]].appTime;

                // update moving average elective appointment waiting time
                if(patient->patientType == 1){
                    if(previousWeek < week[i]){
                        movingAvgElectiveAppWT[previousWeek] = movingAvgElectiveAppWT[previousWeek] / numberOfElectivePerWeek;
                        numberOfElectivePerWeek = 0;
                        previousWeek = week[i];
                    }
                    wt = patient->getAppWT();
                    movingAvgElectiveAppWT[week[i]] += wt;
                    numberOfElectivePerWeek++;
                    avgElectiveAppWT += wt;
                    numberOfElective++;
                }

                // set next slot of the current patient type
                found = false; int startD = day[i]; int startS = slot[i] + 1;
                for(w = week[i]; w < W && !found; w++){
                    for(d = startD; d < D && !found; d++){
                        for(s = startS; s < S && !found; s++){
                            if(weekSchedule[d][s].patientType == patient->patientType){
                                week[i] = w;
                                day[i] = d;
                                slot[i] = s;
                                found = true;
                            }
                        }
                        startS = 0;
                    }
                    startD = 0;
                }
                if(!found) week[i] = W;
            }
        }

        // update moving average elective appointment waiting time in last week
        movingAvgElectiveAppWT[W-1] = movingAvgElectiveAppWT[W-1] / numberOfElectivePerWeek;

        // calculate objective value
        avgElectiveAppWT = avgElectiveAppWT / numberOfElective;
    }

    public void sortPatientsOnAppTime(){
        patients.sort([](const Patient &patient1, const Patient &patient2){
            // unplanned patients at the end of the list in order of their call
            if(patient1.scanWeek == -1 && patient2.scanWeek == -1){
                if (patient1.callWeek != patient2.callWeek)
                    return patient1.callWeek < patient2.callWeek;
                if (patient1.callDay != patient2.callDay)
                    return patient1.callDay < patient2.callDay;
                if (patient1.callTime != patient2.callTime)
                    return patient1.callTime < patient2.callTime;
                if (patient1.scanType == 2)                             // if arrival time same, urgent patient before elective patient
                    return true;
                if (patient2.scanType == 2)
                    return false;
                return true;
            }
            if(patient1.scanWeek == -1){
                return false;
            }
            if(patient2.scanWeek == -1){
                return true;
            }

            if (patient1.scanWeek != patient2.scanWeek)
                return patient1.scanWeek < patient2.scanWeek;
            if (patient1.scanDay != patient2.scanDay)
                return patient1.scanDay < patient2.scanDay;
            if (patient1.appTime != patient2.appTime)
                return patient1.appTime < patient2.appTime;
            if (patient1.scanType == 2)                             // if arrival time same, urgent patient before elective patient
                return true;
            if (patient2.scanType == 2)
                return false;
            if(patient1.nr < patient2.nr){
                return true;
            }
            if(patient1.nr > patient2.nr){
                return false;
            }
            return true;
        });
    }

}
