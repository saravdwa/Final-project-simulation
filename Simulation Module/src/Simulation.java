import java.io.*;
import java.util.*;

import static java.lang.Double.max;

public class Simulation {
    //files
    String inputFileName = "C:\\Users\\info\\IdeaProjects\\Final-project-simulation\\Simulation Module\\src\\input-S1-14.txt";
    String outputFileName = "C:\\Users\\info\\IdeaProjects\\Final-project-simulation\\Simulation Module\\src\\output.txt";


    // Variables and parameters (given in assignment explanation)
    int D = 6;                         // number of days per week (NOTE: Sunday not included! so do NOT use to calculate appointment waiting time)
    int amountOTSlotsPerDay = 10;      // number of overtime slots per day
    int S = 32 + amountOTSlotsPerDay;  // number of slots per day
    double slotLenght = 15.0 / 60.0;             // duration of a slot (in hours)
    double lambdaElective = 28.345; //number of elective arrivals per day follows a Poisson distribution with this lambda
    double meanTardiness = 0;  //tardiness follows a normal distribution with a mean of 0min
    double stdevTardiness = 2.5; //tardiness follows a normal distribution with a standard deviation of 2.5min
    double probNoShow = 0.02; // probability of not showing up is 2%
    double meanElectiveDuration = 15; //planned scans follow a normal distribution with mean of 15min
    double stdevElectiveDuration = 3; //planned scans follow a normal distribution with stdev of 3min
    double[] lambdaUrgent = new double[]{2.5, 1.25}; //urgent patients follow Poisson distribution with lamda = 2.5 (full day) vs. lambda = 1.25 (half a day)
    double[] probUrgentType = new double[]{0.7, 0.1, 0.1, 0.05, 0.05}; //frequencies of scan types urgent patients (see table 2)
    double[] cumulativeProbUrgentType = new double[]{0.7, 0.8, 0.9, 0.95, 1.0};
    double[] meanUrgentDuration = new double[]{15, 17.5, 22.5, 30, 30};
    double[] stdevUrgentDuration = new double[]{2.5, 1, 2.5, 1, 4.5};
    double weightEl = 1.0 / 168.0;      // objective weight elective appointment wait time: 168h/week for planned patients (see objective function)
    double weightUr = 1.0 / 9.0;        // objective weight urgent scan wait time: 9h/week for urgent patients (see objective function)

    int W, R;                          // number of weeks 😊 runs lenght) and number of replications (set their values yourself in the initalization method!)
    int d, s, w, r;
    int rule;                          // the appointment scheduling rule
    Slot[][] weekSchedule;  // array of the cyclic slot schedule (days-slots)


    // Variables within ONE simuation
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

    // Initialization of a "simulation" object
    public Simulation() {
        // Set test case variables
        //TODO: each time you make a different simulation: set these variables to the correct values
            /* we have to make our own input file with the schedule we want to use in the simulation and
               and change the path here (this file is then used in the method runSimulation down below
            */
        W = 10;                      // number of weeks to simulate = run length
        R = 1;                      // number of replications
        rule = 1;                   // the appointment scheduling rule to apply

        // Initialize variables
        avgElectiveAppWT = 0; //in the beginning of the simulation, there are not yet waiting times
        avgElectiveScanWT = 0;
        avgUrgentScanWT = 0;
        avgOT = 0; //no overtime yet
        numberOfElectivePatientsPlanned = 0; //nothing planned, this is only the initialization phase = creating structure
        numberOfUrgentPatientsPlanned = 0;

        // Initialize arrays
        weekSchedule = new Slot[D][s]; //weekschedule has D(=6) rows for the days of the week and s columns with the number of slots per day
        for (this.d = 0; this.d < D; this.d++) {
            weekSchedule[this.d] = new Slot[S]; //each day in the weekschedule is filled with a number of slots
        }
        movingAvgElectiveAppWT = new double[W]; //make an array for the moving average of the waiting time for planned patients with length W
        movingAvgElectiveScanWT = new double[W];
        movingAvgUrgentScanWT = new double[W];
        movingAvgOT = new double[W];
    }

    /*
    This method is called upon in the main method, after creation of a simulation with some parameters
    we have to choose ourselves for each simulation in the constructor (the method before this one)
     */
    public void runSimulations() throws IOException {
        double electiveAppWT = 0;
        double electiveScanWT = 0;
        double urgentScanWT = 0;
        double OT = 0; //overtime
        double OV = 0; //overtime weighted
        setWeekSchedule();
        // set cyclic slot schedule based on given input file
        Random r = new Random();
        System.out.printf("r \t elAppWT \t elScanWT \t urScanWT \t OT \t OV \n");
        // run R replications over all the slots s
        for (s = 0; s < R; s++) {
            resetSystem();          // reset all variables related to 1 replication
            r.setSeed(s); // set seed value for random value generator
            runOneSimulation();     // run 1 simulation / replication
            electiveAppWT += avgElectiveAppWT;
            electiveScanWT += avgElectiveScanWT;
            urgentScanWT += avgUrgentScanWT;
            OT += avgOT;
            OV += avgElectiveAppWT / weightEl + (avgUrgentScanWT / weightUr);
            System.out.printf("%d \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \n", r, avgElectiveAppWT, avgElectiveScanWT, avgUrgentScanWT, avgOT, avgElectiveAppWT / weightEl + avgUrgentScanWT / weightUr);
        }
        electiveAppWT = electiveAppWT / R; //divide by the number of replications to get the value once
        electiveScanWT = electiveScanWT / R;
        urgentScanWT = urgentScanWT / R;
        OT = OT / R;
        OV = OV / R;
        double objectiveValue = electiveAppWT / weightEl + urgentScanWT / weightUr;
        System.out.printf("Avg.: \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \n", electiveAppWT, electiveScanWT, urgentScanWT, OT, objectiveValue);

        // print results
        // FILE *file = fopen("/Users/tinemeersman/Documents/project SMA 2022 student code /output.txt", "a"); // TODO: use your own directory
        // TODO: print the output you need to a .txt file: @Emma kan jij dit doen?
        //fclose(file);


        File file = new File(outputFileName);               // if file doesnt exists, then create it
        if (!file.exists()) {
            file.createNewFile();                           // create the file
        } else {
            PrintWriter writer = new PrintWriter(file);     // empty the file
            writer.print("");
            writer.close();
        }
        FileWriter fileWriter = new FileWriter(file.getAbsoluteFile(), true); // APPENDS the text file with anything printed to the file during the rest of the procedure
        PrintWriter printWriter = new PrintWriter(fileWriter);                       // OPEN OUTPUT FILE
        printWriter.println("test");
    }

    public void setWeekSchedule() {
        // Read the schedule from the textfile and set the slot types (0=none, 1=elective, 2=urgent within normal working hours)
        InputStream input = null;
        try {
            input = new FileInputStream(inputFileName);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
        Scanner inputFile = new Scanner(input);
        int elementInt;
        while (inputFile.hasNext()) {
            for (s = 0; s < 32; s++) {
                for (d = 0; d < D; d++) {
                    elementInt = inputFile.nextInt();
                    weekSchedule[d][s].slotType = elementInt;
                    weekSchedule[d][s].patientType = elementInt;
                }
            }
        }

        inputFile.close();

        // Set the type of the overtime slots (3=urgent in overtime)
        for (d = 0; d < D; d++) {
            for (s = 32; s < S; s++) {
                weekSchedule[d][s].slotType = 3;
                weekSchedule[d][s].patientType = 2;
            }
        }

        // set start and appointment time
        double time;
        for (d = 0; d < D; d++) {
            time = 8; // start time slot schedule
            int countPatient = 0;
            int B = 2; // 2 blocks for blocking rule
            for (s = 0; s < S; s++) {
                // start time slot
                weekSchedule[d][s].startTime = time;
                // appointment time slot
                if (weekSchedule[d][s].slotType != 1) {    // all slot types not elective : appointment time = slot start time (because urgent!)
                    weekSchedule[d][s].appTime = time;
                } else {                                   // elective slots: appointment time is set according to rule ! (because can be late)
                    if (rule == 1) { // FIFO
                        weekSchedule[d][s].appTime = time;
                    } else if (rule == 2) {
                        // TODO: Bailey-Welch rule
                        if (countPatient < 2) {
                            weekSchedule[d][countPatient].appTime = time; //should the first two patients have a different time slot, both starting at 8 or the same one?
                        } else {
                            weekSchedule[d][s].appTime = time - slotLenght;
                        }
                        countPatient++;
                    } else if (rule == 3) {
                        // TODO: Blocking rule
                        if (s < (S / B)) {  // should this be programmed with the idea in mind that B can be other than 2?
                            weekSchedule[d][s].appTime = 8; //first block always starts at 8
                        } else {
                            weekSchedule[d][s].appTime = 8 + slotLenght*(S/B);
                        }
                    } else if (rule == 4) {
                        // TODO: Benchmark rule
                        double k = 0.5;
                        weekSchedule[d][s].appTime = time - k*stdevElectiveDuration;
                    }
                }

                //update time variable
                time += slotLenght;
                if (time == 12) {
                    time = 13;
                } // skip to the end of the lunch break
            }
        }

    }

    public void resetSystem() {
        // reset all variables related to 1 replication
        patients.clear();
        avgElectiveAppWT = 0;
        avgElectiveScanWT = 0;
        avgUrgentScanWT = 0;
        avgOT = 0;
        numberOfElectivePatientsPlanned = 0;
        numberOfUrgentPatientsPlanned = 0;

        for (w = 0; w < W; w++) {
            movingAvgElectiveAppWT[w] = 0;
            movingAvgElectiveScanWT[w] = 0;
            movingAvgUrgentScanWT[w] = 0;
            movingAvgOT[w] = 0;
        }
    }

    void runOneSimulation() {
        generatePatients();     // create patient arrival events (elective patients call, urgent patient arrive at the hospital)
        schedulePatients();     // schedule urgent and elective patients in slots based on their arrival events => detrmine the appointment wait time
        sortPatientsOnAppTime();   // sort patients on their appointment time (unscheduled patients are grouped at the end of the list)

        // determine scan wait time per patient and overtime per day
        int prevWeek = 0;
        int prevDay = -1;
        int[] numberOfPatientsWeek = new int[]{0, 0};
        int[] numberOfPatients = new int[]{0, 0};
        double arrivalTime, wt;
        double prevScanEndTime = 0;
        boolean prevIsNoShow = false;
        // go over arrival events (i.e. the moment the patient arrives at the hospital)
        for (int i = 0; i < patients.size(); i++) {

            if (patients.get(i).scanWeek == -1) { // stop at the first unplanned patient
                break;
            }


            arrivalTime = (double) patients.get(i).appTime + patients.get(i).tardiness;
            // SCAN WT
            if (!patients.get(i).isNoShow) {
                if (patients.get(i).scanWeek != prevWeek || patients.get(i).scanDay != prevDay) {
                    patients.get(i).scanTime = arrivalTime;
                } else {
                    if (prevIsNoShow) {
                        patients.get(i).scanTime = max(weekSchedule[patients.get(i).scanDay][patients.get(i).slotNr].startTime, max(prevScanEndTime, arrivalTime)); // note we assume we wait at least 15minutes on a no-show patient to see whether he shows or is just late
                    } else {
                        patients.get(i).scanTime = max(prevScanEndTime, arrivalTime);
                    }
                }
                wt = patients.get(i).getScanWT();
                if (patients.get(i).patientType == 1) {
                    movingAvgElectiveScanWT[patients.get(i).scanWeek] += wt;
                } else {
                    movingAvgUrgentScanWT[patients.get(i).scanWeek] += wt;
                }
                numberOfPatientsWeek[patients.get(i).patientType - 1]++;
                if (patients.get(i).patientType == 1) {
                    avgElectiveScanWT += wt;
                } else {
                    avgUrgentScanWT += wt;
                }
                numberOfPatients[patients.get(i).patientType - 1]++;
            }

            // OVERTIME
            if (prevDay > -1 && prevDay != patients.get(i).scanDay) {
                if (d == 3 || d == 5) {
                    movingAvgOT[prevWeek] += max(0.0, prevScanEndTime - 13);
                } else {
                    movingAvgOT[prevWeek] += max(0.0, prevScanEndTime - 17);
                }
                if (d == 3 || d == 5) {
                    avgOT += max(0.0, prevScanEndTime - 13);
                } else {
                    avgOT += max(0.0, prevScanEndTime - 17);
                }
            }

            // update moving averages if week ends
            if (prevWeek != patients.get(i).scanWeek) {
                movingAvgElectiveScanWT[prevWeek] = movingAvgElectiveScanWT[prevWeek] / numberOfPatientsWeek[0];
                movingAvgUrgentScanWT[prevWeek] = movingAvgUrgentScanWT[prevWeek] / numberOfPatientsWeek[1];
                movingAvgOT[prevWeek] = movingAvgOT[prevWeek] / D;
                numberOfPatientsWeek[0] = 0;
                numberOfPatientsWeek[1] = 0;
            }

            //set prev patient
            if (patients.get(i).isNoShow) {
                //prevScanEndTime stays the same, it is the end time of the patient before the no-show patient
                prevIsNoShow = true;
            } else {
                prevScanEndTime = patients.get(i).scanTime + patients.get(i).duration;
                prevIsNoShow = false;
            }
            prevWeek = patients.get(i).scanWeek;
            prevDay = patients.get(i).scanDay;
        }
        // update moving averages of the last week
        movingAvgElectiveScanWT[W - 1] = movingAvgElectiveScanWT[W - 1] / numberOfPatientsWeek[0];
        movingAvgUrgentScanWT[W - 1] = movingAvgUrgentScanWT[W - 1] / numberOfPatientsWeek[1];
        movingAvgOT[W - 1] = movingAvgOT[W - 1] / D;

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


    public int getRandomScanType() {

        Random s = new Random();

        float r;
        r = (float) s.nextInt(1000 - 1) / 1000;

        int type = -1;
        for (int i = 0; i < 5 && type == -1; i++) {
            if (r < cumulativeProbUrgentType[i]) {
                type = i;
            }
        }
        return type;
    }

    public void generatePatients() {
        double arrivalTimeNext;
        int counter = 0; // total number of patients so far
        int patientType, scanType, endTime;
        double callTime, tardiness, duration, lambda;
        boolean noShow;
        Random r = null;
        for (w = 0; w < W; w++) {
            for (d = 0; d < D; d++) { // not on Sunday
                // generate ELECTIVE patients for this day
                if (d < D - 1) {  // not on Saturday either
                    arrivalTimeNext = 8 + Exponential_distribution(lambdaElective, r) * (17 - 8);
                    while (arrivalTimeNext < 17) { // desk open from 8h until 17h
                        patientType = 1;                // elective
                        scanType = 0;                   // no scan type
                        callTime = arrivalTimeNext;     // set call time, i.e. arrival event time
                        tardiness = Normal_distribution(meanTardiness, stdevTardiness, r) / 60.0;       // in practice this is not known yet at time of call
                        noShow = Bernouilli_distribution(probNoShow, r);                                // in practice this is not known yet at time of call
                        duration = Normal_distribution(meanElectiveDuration, stdevElectiveDuration, r) / 60.0; // in practice this is not known yet at time of call
                        Patient patient = new Patient(counter, patientType, scanType, w, d, callTime, tardiness, noShow, duration);
                        patients.add(patient);
                        counter++;
                        arrivalTimeNext = arrivalTimeNext + Exponential_distribution(lambdaElective, r) * (17 - 8); // arrival time of next patient (if < 17h)
                    }
                }

                // generate URGENT patients for this day
                if (d == 3 || d == 5) {
                    lambda = lambdaUrgent[1]; // on Wed and Sat, only half a day!
                    endTime = 12;
                } else {
                    lambda = lambdaUrgent[0];
                    endTime = 17;
                }
                arrivalTimeNext = 8 + Exponential_distribution(lambda, r) * (endTime - 8);
                while (arrivalTimeNext < endTime) { // desk open from 8h until 17h
                    patientType = 2;                // urgent
                    scanType = getRandomScanType(); // set scan type
                    callTime = arrivalTimeNext;     // set arrival time, i.e. arrival event time
                    tardiness = 0;                  // urgent patients have an arrival time = arrival event time
                    noShow = false;                 // urgent patients are never no-show
                    duration = Normal_distribution(meanUrgentDuration[scanType], stdevUrgentDuration[scanType], r) / 60.0; // in practice this is not known yet at time of arrival
                    Patient patient = new Patient(counter, patientType, scanType, w, d, callTime, tardiness, noShow, duration);
                    //ArrayList<patient> pa = new ArrayList<patient>();
                    patients.add(patient);
                    counter++;
                    arrivalTimeNext = arrivalTimeNext + Exponential_distribution(lambda, r) * (endTime - 8); // arrival time of next patient (if < 17h)
                }
            }
        }
    }

    private boolean Bernouilli_distribution(double probNoShow, Random r) {
        double j1 = (float) r.nextInt(1000 + 1) / 1000;
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
        do {
            v1 = (float) r.nextInt(1000 + 1) * 2;
            v1 /= 1000;
            v1 -= 1;
            v2 = (float) r.nextInt(1000 + 1) * 2;
            v2 /= 1000;
            v2 -= 1;
            t = v1 * v1 + v2 * v2;
        }
        while (t >= 1 || t == 0);
        double multiplier = Math.sqrt(-2 * Math.log(t) / t);
        return (int) (v1 * multiplier * stdevTardiness + meanTardiness);

    }

    private double Exponential_distribution(double lambdaElective, Random r) {
        double j1 = (float) r.nextInt(1000 + 1) / 1000;
        if (j1 == 0)
            j1 += 0.0001;
        double j2 = -Math.log(j1) / lambdaElective;

        return j2;
    }

    public int getNextSlotNrFromTime(int day, int patientType, double time) {
        boolean found = false;
        int slotNr = -1;
        for (s = 0; !found && s < S; s++) {
            if (weekSchedule[day][s].appTime > time && patientType == weekSchedule[day][s].patientType) {
                found = true;
                slotNr = s;
            }
        }
        if (!found) {
            System.out.printf("NO SLOT EXISTS DURING TIME %.2f \n", time);
            System.exit(0);
        }
        return slotNr;
    }

    public void schedulePatients1() {
        Collections.sort(patients, new Comparator<Patient>() {

            public int compare(Patient p1, Patient p2) {
                Integer weekCount1 = p1.callWeek;
                Integer weekCount2 = p2.callWeek;
                Integer dayCount1 = p1.callDay;
                Integer dayCount2 = p2.callDay;
                Double timeCount1 = p1.callTime;
                Double timeCount2 = p2.callTime;

                if (weekCount1.compareTo(weekCount2) != 0) {
                    return (weekCount1 < weekCount2) ? 1 : 0; //1 if true, anders 0
                }
                if (dayCount1.compareTo(dayCount2) != 0) {
                    return (dayCount1 < dayCount2) ? 1 : 0;
                }
                if (timeCount1.compareTo(timeCount2) != 0) {
                    return (timeCount1 < timeCount2) ? 1 : 0;
                }
                if (p1.scanType == 2)
                    return 1;
                if (p2.scanType == 2)
                    return 0;

                return 1;
            }
        });

        int[] week = {0, 0}; // week of the next available slot {elective,urgent}
        int[] day = {0, 0}; // day of the next available slot {elective,urgent}
        int[] slot = {0, 0}; // slotNr of the next available slot {elective,urgent}

        //find first slot of each patient type (note, we assume each day (i.e. also day 0) has at least one slot of each patient type!)
        //elective
        d = 0;
        boolean found = false;
        for (s = 0; s < S && !found; s++) {
            if (weekSchedule[d][s].patientType == 1) {
                day[0] = d;
                slot[0] = s;
                found = true;
            }
            //urgent
            found = false;
            for (s = 0; s < S && !found; s++) {
                if (weekSchedule[d][s].patientType == 2) {
                    day[1] = d;
                    slot[1] = s;
                    found = true;
                }
            }
        }

    }

    public void schedulePatients() {
        Collections.sort(patients, new Comparator<Patient>() {

            public int compare(Patient p1, Patient p2) {
                Integer weekCount1 = p1.callWeek;
                Integer weekCount2 = p2.callWeek;
                Integer dayCount1 = p1.callDay;
                Integer dayCount2 = p2.callDay;
                Double timeCount1 = p1.callTime;
                Double timeCount2 = p2.callTime;

                if (weekCount1.compareTo(weekCount2) != 0) {
                    return (weekCount1 < weekCount2) ? 1 : 0; //1 if true, anders 0
                }
                if (dayCount1.compareTo(dayCount2) != 0) {
                    return (dayCount1 < dayCount2) ? 1 : 0;
                }
                if (timeCount1.compareTo(timeCount2) != 0) {
                    return (timeCount1 < timeCount2) ? 1 : 0;
                }
                if (p1.scanType == 2)
                    return 1;
                if (p2.scanType == 2)
                    return 0;

                return 1;
            }
        });

        int[] week = {0, 0}; // week of the next available slot {elective,urgent}
        int[] day = {0, 0}; // day of the next available slot {elective,urgent}
        int[] slot = {0, 0}; // slotNr of the next available slot {elective,urgent}

        //find first slot of each patient type (note, we assume each day (i.e. also day 0) has at least one slot of each patient type!)
        //elective
        d = 0;
        boolean found = false;
        for (s = 0; s < S && !found; s++) {
            if (weekSchedule[d][s].patientType == 1) {
                day[0] = d;
                slot[0] = s;
                found = true;
            }
            //urgent
            found = false;
            for (s = 0; s < S && !found; s++) {
                if (weekSchedule[d][s].patientType == 2) {
                    day[1] = d;
                    slot[1] = s;
                    found = true;
                }
            }
        }

        // go over SORTED patient list and assign slots
        int previousWeek = 0;
        int numberOfElective = 0;
        int numberOfElectivePerWeek = 0;   // keep track of week to know when to update moving average elective appointment waiting time
        double wt;
        int slotNr;
        for (int i = 0; i < patients.size(); i++) {
            //Patient *pat = &*patient;

            //set index i dependant on patient type
            i = patients.get(i).patientType - 1;

            // if still within the planning horizon:
            if (week[i] < W) {

                // determine week where we start searching for a slot
                if (patients.get(i).callWeek > week[i]) {
                    week[i] = patients.get(i).callWeek;
                    day[i] = 0;
                    slot[i] = getNextSlotNrFromTime(day[i], patients.get(i).patientType, 0);           // note we assume there is at least one slot of each patient type per day => this line will find first slot of this type
                }
                // determine day where we start searching for a slot
                if (patients.get(i).callWeek == week[i] && patients.get(i).callDay > day[i]) {
                    day[i] = patients.get(i).callDay;
                    slot[i] = getNextSlotNrFromTime(day[i], patients.get(i).patientType, 0);           // note we assume there is at least one slot of each patient type per day => this line will find first slot of this type
                }
                // determine slot
                if (patients.get(i).callWeek == week[i] && patients.get(i).callDay == day[i] && patients.get(i).callTime >= weekSchedule[day[i]][slot[i]].appTime) {
                    // find last slot on day "day[i]"
                    found = false;
                    slotNr = -1;
                    for (s = S - 1; s >= 0 && !found; s--) {
                        if (weekSchedule[day[i]][s].patientType == patients.get(i).patientType) {
                            found = true;
                            slotNr = s;
                        }
                    }
                    // urgent patients have to be treated on the same day either in normal hours or in overtime (!! make sure there are enough overtime slots)
                    // for elective patients: check if the patient call time is before the last slot, i.e. if the patient can be planned on day "day[i]"
                    if (patients.get(i).patientType == 2 || patients.get(i).callTime < weekSchedule[day[i]][slotNr].appTime) {
                        slot[i] = getNextSlotNrFromTime(day[i], patients.get(i).patientType, patients.get(i).callTime);   // find the first elective slot after the call time on day "day[i]"
                    } else {
                        // determine the next day
                        if (day[i] < D - 1) {
                            day[i] = day[i] + 1;
                        } else {
                            day[i] = 0;
                            week[i] = week[i] + 1;
                        }
                        if (week[i] < W) {   // find the first slot on the next day (if within the planning horizon)
                            slot[i] = getNextSlotNrFromTime(day[i], patients.get(i).patientType, 0);
                        }
                    }
                }

                // schedule the patient
                patients.get(i).scanWeek = week[i];
                patients.get(i).scanDay = day[i];
                patients.get(i).slotNr = slot[i];
                patients.get(i).appTime = weekSchedule[day[i]][slot[i]].appTime;

                // update moving average elective appointment waiting time
                if (patients.get(i).patientType == 1) {
                    if (previousWeek < week[i]) {
                        movingAvgElectiveAppWT[previousWeek] = movingAvgElectiveAppWT[previousWeek] / numberOfElectivePerWeek;
                        numberOfElectivePerWeek = 0;
                        previousWeek = week[i];
                    }
                    wt = patients.get(i).getAppWT();
                    movingAvgElectiveAppWT[week[i]] += wt;
                    numberOfElectivePerWeek++;
                    avgElectiveAppWT += wt;
                    numberOfElective++;
                }

                // set next slot of the current patient type
                found = false;
                int startD = day[i];
                int startS = slot[i] + 1;
                for (w = week[i]; w < W && !found; w++) {
                    for (d = startD; d < D && !found; d++) {
                        for (s = startS; s < S && !found; s++) {
                            if (weekSchedule[d][s].patientType == patients.get(i).patientType) {
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
                if (!found) week[i] = W;
            }
        }

        // update moving average elective appointment waiting time in last week
        movingAvgElectiveAppWT[W - 1] = movingAvgElectiveAppWT[W - 1] / numberOfElectivePerWeek;

        // calculate objective value
        avgElectiveAppWT = avgElectiveAppWT / numberOfElective;
    }

    public void sortPatientsOnAppTime() {
        Collections.sort(patients, (p1, p2) -> {
            Integer weekCount1 = p1.callWeek;
            Integer weekCount2 = p2.callWeek;
            Integer dayCount1 = p1.callDay;
            Integer dayCount2 = p2.callDay;
            Double timeCount1 = p1.callTime;
            Double timeCount2 = p2.callTime;
            Integer scanWeekCount1 = p1.scanWeek;
            Integer scanWeekCount2 = p2.scanWeek;
            Integer scanDayCount1 = p1.scanDay;
            Integer scanDayCount2 = p2.scanDay;
            ;
            Double scanAppTimeCount1 = p1.appTime;
            Double scanAppTimeCount2 = p2.appTime;

            if (p1.scanWeek == -1 && p2.scanWeek == -1) {
                if (weekCount1.compareTo(weekCount2) != 0) {
                    return (weekCount1 < weekCount2) ? 1 : 0; //1 if true, anders 0
                }
                if (dayCount1.compareTo(dayCount2) != 0) {
                    return (dayCount1 < dayCount2) ? 1 : 0;
                }
                if (timeCount1.compareTo(timeCount2) != 0) {
                    return (timeCount1 < timeCount2) ? 1 : 0;
                }
                if (p1.scanType == 2)
                    return 1;
                if (p2.scanType == 2)
                    return 0;
                return 1;
            }
            if (p1.scanWeek == -1)
                return 0;
            if (p2.scanWeek == -1)
                return 1;
            //return 1: should we exit the function here?
            if (scanWeekCount1.compareTo(scanWeekCount2) != 0) {
                return (scanWeekCount1 < scanWeekCount2) ? 1 : 0;   //1 if true, anders 0
            }
            if (scanDayCount1.compareTo(scanDayCount2) != 0) {
                return (scanDayCount1 < scanDayCount2) ? 1 : 0;
            }
            if (scanAppTimeCount1.compareTo(scanAppTimeCount2) != 0) {
                return (scanAppTimeCount1 < scanAppTimeCount2) ? 1 : 0;
            }
            if (p1.scanType == 2)
                return 1;
            if (p2.scanType == 2)
                return 0;
            if (p1.nr < p2.nr)
                return 1;
            if (p1.nr > p2.nr)
                return 0;
            return 1;
        });
    }


}
