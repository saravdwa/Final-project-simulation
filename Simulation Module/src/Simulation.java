import java.io.*;
import java.util.*;

import static java.lang.Double.max;
import static java.lang.StrictMath.log;
import static java.lang.StrictMath.sqrt;

public class Simulation {


    //files : ik heb een relatief pad gemaakt, zodat het op iedereen zijn computer werkt
    String filePath = "Simulation Module/src/input-S1-14.txt";
    File file = new File(filePath);
    String inputFileName = file.getPath();


    // Variables and parameters (given in assignment explanation)
    int D = 6;                         // number of days per week (NOTE: Sunday not included! so do NOT use to calculate appointment waiting time)
    int amountOTSlotsPerDay = 10;      // number of overtime slots per day
    int S = 32 + amountOTSlotsPerDay;  // number of slots per day
    int countIterationForWelch = 0;
    int warmUpRemoved = 325;
    double slotLength = 15.0 / 60.0;             // duration of a slot (in hours)
    double lambdaElective = 28.345; //number of elective arrivals per day follows a Poisson distribution with this lambda
    double meanTardiness = 0;  //tardiness follows a normal distribution with a mean of 0min
    double stdevTardiness = 2.5; //tardiness follows a normal distribution with a standard deviation of 2.5min
    double probNoShow = 0.02; // probability of not showing up is 2%
    double meanElectiveDuration = 15; //planned scans follow a normal distribution with mean of 15min
    double stdevElectiveDuration = 3; //planned scans follow a normal distribution with stdev of 3min
    double[] lambdaUrgent = new double[]{2.5, 1.25}; //urgent patients follow Poisson distribution with lamda = 2.5 (full day) vs. lambda = 1.25 (half a day)
    //double[] probUrgentType = new double[]{0.7, 0.1, 0.1, 0.05, 0.05}; //frequencies of scan types urgent patients (see table 2)
    double[] cumulativeProbUrgentType = new double[]{0.7, 0.8, 0.9, 0.95, 1.0};
    double[] meanUrgentDuration = new double[]{15, 17.5, 22.5, 30, 30};
    double[] stdevUrgentDuration = new double[]{2.5, 1, 2.5, 1, 4.5};
    double weightEl = 1.0 / 168.0;      // objective weight elective appointment wait time: 168h/week for planned patients (see objective function)
    double weightUr = 1.0 / 9.0;        // objective weight urgent scan wait time: 9h/week for urgent patients (see objective function)

    int W , R;                          // number of weeks 😊 runs length) and number of replications (set their values yourself in the initalization method!)
    int d, s, w;
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

    Random random = new Random();
    Random commonR = new Random();
    Random randFct = new Random();
    Random[] listRandom = new Random[451];
    Random[] listAntiRandom = new Random[451];

    public Random generateRandom(){
        return randFct;
    }


    // Initialization of a "simulation" object
    public Simulation() {
        // Set test case variables
        //TODO: each time you make a different simulation: set these variables to the correct values
        W = 100;                      // number of weeks to simulate = run length
        R = 450;                      // number of replications
        rule = 3;                   // the appointment scheduling rule to apply

        // Initialize variables
        avgElectiveAppWT = 0; //in the beginning of the simulation, there are not yet waiting times
        avgElectiveScanWT = 0;
        avgUrgentScanWT = 0;
        avgOT = 0; //no overtime yet
        numberOfElectivePatientsPlanned = 0; //nothing planned, this is only the initialization phase = creating structure
        numberOfUrgentPatientsPlanned = 0;

        // Initialize arrays
        weekSchedule = new Slot[D][S]; //weekschedule has D(=6) rows for the days of the week and s columns with the number of slots per day
        for (d = 0; d < D; d++) {
            for (int s = 0; s < 32; s++) {
                weekSchedule[d][s] = new Slot();
            }
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
        double OV2 = 0;
        setWeekSchedule();
        // set cyclic slot schedule based on given input file
        System.out.print("r \t common random number \t elAppWT \t elScanWT \t urScanWT \t OT \t \t \t OV \n ");

        // run R replications over all the slots s
        for (int f = 1; f <= R; f++) {
            resetSystem();          // reset all variables related to 1 replication
            randFct.setSeed(f);
            listRandom[f] = randFct;
            //listAntiRandom[f].nextFloat() = 1 - randFct.nextFloat();
            commonR.setSeed(f);// set seed value for random value generator
            runOneSimulation();     // run 1 simulation / replication
            countIterationForWelch++;
            if (countIterationForWelch > 125) {
                electiveAppWT += avgElectiveAppWT;
                electiveScanWT += avgElectiveScanWT;
                urgentScanWT += avgUrgentScanWT;
                OT += avgOT;
                OV += avgElectiveAppWT * weightEl + (avgUrgentScanWT * weightUr);
                OV2 += 5;
                float getRandom = commonR.nextFloat();
                System.out.printf("%d \t %.6f \t \t\t %.6f \t %.6f \t \t %.6f \t %.6f \t %.6f \n", f, getRandom, avgElectiveAppWT, avgElectiveScanWT, avgUrgentScanWT, avgOT, avgElectiveAppWT * weightEl + avgUrgentScanWT * weightUr);

            }
        }

        electiveAppWT = electiveAppWT / warmUpRemoved; //divide by the number of replications to get the value once
        electiveScanWT = electiveScanWT / warmUpRemoved;
        urgentScanWT = urgentScanWT / warmUpRemoved;
        OT = OT / warmUpRemoved;
        OV = OV / warmUpRemoved;
        double objectiveValue = electiveAppWT * weightEl + urgentScanWT * weightUr;
        System.out.printf("Avg.: \t \t \t \t\t\t %.6f \t %.6f \t\t %.6f \t %.6f \t %.6f \n", electiveAppWT, electiveScanWT, urgentScanWT, OT, objectiveValue);
    }

    public void setWeekSchedule() throws FileNotFoundException {
        // Read the schedule from the text-file and set the slot types (0=none, 1=elective, 2=urgent within normal working hours)
        InputStream input = new FileInputStream(inputFileName);
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
                weekSchedule[d][s] = new Slot();
                weekSchedule[d][s].slotType = 3;
                weekSchedule[d][s].patientType = 2;
            }
        }

        // set start and appointment time
        double time;
        double startDay = 8;
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
                    if (rule == 1) {
                        // FIFO
                        weekSchedule[d][s].appTime = time;
                    } else if (rule == 2) {
                        //  Bailey-Welch rule
                        if (countPatient < 2) {
                            weekSchedule[d][countPatient].appTime = startDay;
                        } else {
                            weekSchedule[d][s].appTime = time - slotLength;
                        }
                        countPatient++;
                    } else if (rule == 3) {
                        // Blocking rule
                        if (s % 2 == 0) {
                            weekSchedule[d][s].appTime = time;
                        } else {
                            weekSchedule[d][s].appTime = time - slotLength;
                        }

                    } else if (rule == 4) {
                        // Benchmark rule
                        double k = 0.5;
                        weekSchedule[d][s].appTime = time - k * stdevElectiveDuration;
                    }
                }

                //update time variable
                time += slotLength;
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
        for (Patient patient : patients) {

            if (patient.scanWeek == -1) { // stop at the first unplanned patient
                break;
            }


            arrivalTime = patient.appTime + patient.tardiness;
            // SCAN WT
            if (!patient.isNoShow) {
                if (patient.scanWeek != prevWeek || patient.scanDay != prevDay) {
                    patient.scanTime = arrivalTime;
                } else {
                    if (prevIsNoShow) {
                        patient.scanTime = max(weekSchedule[patient.scanDay][patient.slotNr].startTime, max(prevScanEndTime, arrivalTime)); // note we assume we wait at least 15minutes on a no-show patient to see whether he shows or is just late
                    } else {
                        patient.scanTime = max(prevScanEndTime, arrivalTime);
                    }
                }
                wt = patient.getScanWT();
                if (patient.patientType == 1) {
                    movingAvgElectiveScanWT[patient.scanWeek] += wt;
                } else {
                    movingAvgUrgentScanWT[patient.scanWeek] += wt;
                }
                numberOfPatientsWeek[patient.patientType - 1]++;
                if (patient.patientType == 1) {
                    avgElectiveScanWT += wt;
                } else {
                    avgUrgentScanWT += wt;
                }
                numberOfPatients[patient.patientType - 1]++;
            }

            // OVERTIME
            if (prevDay > -1 && prevDay != patient.scanDay) {
                if (prevDay == 3 || prevDay == 5) {
                    movingAvgOT[prevWeek] += max(0.0, prevScanEndTime - 13);
                } else {
                    movingAvgOT[prevWeek] += max(0.0, prevScanEndTime - 17);
                }
                if (prevDay == 3 || prevDay == 5) {
                    avgOT += max(0.0, prevScanEndTime - 13);
                } else {
                    avgOT += max(0.0, prevScanEndTime - 17);
                }
            }

            // update moving averages if week ends
            if (prevWeek != patient.scanWeek) {
                movingAvgElectiveScanWT[prevWeek] = movingAvgElectiveScanWT[prevWeek] / numberOfPatientsWeek[0];
                movingAvgUrgentScanWT[prevWeek] = movingAvgUrgentScanWT[prevWeek] / numberOfPatientsWeek[1];
                movingAvgOT[prevWeek] = movingAvgOT[prevWeek] / D;
                numberOfPatientsWeek[0] = 0;
                numberOfPatientsWeek[1] = 0;
            }

            //set prev patient
            if (patient.isNoShow) {
                prevIsNoShow = true;
                if (patient.scanWeek != prevWeek || patient.scanDay != prevDay) {

                    prevScanEndTime = weekSchedule[patient.scanDay][patient.slotNr].startTime;

                }
                //prevScanEndTime stays the same, it is the end time of the patient before the no-show patient

            } else {
                prevScanEndTime = patient.scanTime + patient.duration;
                prevIsNoShow = false;
            }
            prevWeek = patient.scanWeek;
            prevDay = patient.scanDay;
        }
        // update moving averages of the last week
        movingAvgElectiveScanWT[W - 1] = movingAvgElectiveScanWT[W - 1] / numberOfPatientsWeek[0];
        movingAvgUrgentScanWT[W - 1] = movingAvgUrgentScanWT[W - 1] / numberOfPatientsWeek[1];
        movingAvgOT[W - 1] = movingAvgOT[W - 1] / D;

        // calculate objective values
        avgElectiveScanWT = avgElectiveScanWT / numberOfPatients[0];
        avgUrgentScanWT = avgUrgentScanWT / numberOfPatients[1];
        avgOT = avgOT / (D * W);

    }


    public int getRandomScanType() {
        float r;
        r = (float) random.nextInt(1000 - 1) / 1000;

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
        for (w = 0; w < W; w++) {
            for (d = 0; d < D; d++) { // not on Sunday
                // generate ELECTIVE patients for this day
                if (d < D - 1) {  // not on Saturday either
                    arrivalTimeNext = 8 + exponential_distribution(lambdaElective) * (17 - 8);
                    while (arrivalTimeNext < 17) { // desk open from 8h until 17h

                        patientType = 1;                // elective
                        scanType = 0;                   // no scan type
                        callTime = arrivalTimeNext;     // set call time, i.e. arrival event time
                        tardiness = normal_distribution(meanTardiness, stdevTardiness) / 60.0;       // in practice this is not known yet at time of call

                        noShow = bernouilli_distribution(probNoShow);                                // in practice this is not known yet at time of call
                        duration = normal_distribution(meanElectiveDuration, stdevElectiveDuration) / 60.0; // in practice this is not known yet at time of call
                        Patient patient = new Patient(counter, patientType, scanType, w, d, callTime, tardiness, noShow, duration);
                        patients.add(patient);
                        counter++;
                        arrivalTimeNext = arrivalTimeNext + exponential_distribution(lambdaElective) * (17 - 8); // arrival time of next patient (if < 17h)
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
                arrivalTimeNext = 8 + exponential_distribution(lambda) * (endTime - 8);
                while (arrivalTimeNext < endTime) { // desk open from 8h until 17h
                    patientType = 2;                // urgent
                    scanType = getRandomScanType(); // set scan type
                    callTime = arrivalTimeNext;     // set arrival time, i.e. arrival event time
                    tardiness = 0;                  // urgent patients have an arrival time = arrival event time
                    noShow = false;                 // urgent patients are never no-show
                    duration = normal_distribution(meanUrgentDuration[scanType], stdevUrgentDuration[scanType]) / 60.0; // in practice this is not known yet at time of arrival
                    Patient patient = new Patient(counter, patientType, scanType, w, d, callTime, tardiness, noShow, duration);
                    //ArrayList<patient> pa = new ArrayList<patient>();
                    patients.add(patient);
                    counter++;
                    arrivalTimeNext = arrivalTimeNext + exponential_distribution(lambda) * (endTime - 8); // arrival time of next patient (if < 17h)
                }
            }
        }
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


    public void schedulePatients() {
        patients.sort((p1, p2) -> {
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
            if (p1.patientType == 2)
                return 1;
            if (p2.patientType == 2)
                return 0;

            return 1;
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

        // go over SORTED patient list and assign slots
        int previousWeek = 0;
        int numberOfElective = 0;
        int numberOfElectivePerWeek = 0;   // keep track of week to know when to update moving average elective appointment waiting time
        double wt;
        int slotNr;
        for (
                Patient value : patients) {
            //Patient *pat = &*patient;

            //set index i dependant on patient type
            int i = value.patientType - 1;

            // if still within the planning horizon:
            if (week[i] < W) {

                // determine week where we start searching for a slot
                if (value.callWeek > week[i]) {
                    week[i] = value.callWeek;
                    day[i] = 0;
                    slot[i] = getNextSlotNrFromTime(day[i], value.patientType, 0);           // note we assume there is at least one slot of each patient type per day => this line will find first slot of this type
                }
                // determine day where we start searching for a slot
                if (value.callWeek == week[i] && value.callDay > day[i]) {
                    day[i] = value.callDay;
                    slot[i] = getNextSlotNrFromTime(day[i], value.patientType, 0);           // note we assume there is at least one slot of each patient type per day => this line will find first slot of this type
                }
                // determine slot
                if (value.callWeek == week[i] && value.callDay == day[i] && value.callTime >= weekSchedule[day[i]][slot[i]].appTime) {
                    // find last slot on day "day[i]"
                    found = false;
                    slotNr = -1;
                    for (s = S - 1; s >= 0 && !found; s--) {
                        if (weekSchedule[day[i]][s].patientType == value.patientType) {
                            found = true;
                            slotNr = s;
                        }
                    }
                    // urgent patients have to be treated on the same day either in normal hours or in overtime (!! make sure there are enough overtime slots)
                    // for elective patients: check if the patient call time is before the last slot, i.e. if the patient can be planned on day "day[i]"
                    if (value.patientType == 2 || value.callTime < weekSchedule[day[i]][slotNr].appTime) {
                        slot[i] = getNextSlotNrFromTime(day[i], value.patientType, value.callTime);   // find the first elective slot after the call time on day "day[i]"
                    } else {
                        // determine the next day
                        if (day[i] < D - 1) {
                            day[i] = day[i] + 1;
                        } else {
                            day[i] = 0;
                            week[i] = week[i] + 1;
                        }
                        if (week[i] < W) {   // find the first slot on the next day (if within the planning horizon)
                            slot[i] = getNextSlotNrFromTime(day[i], value.patientType, 0);
                        }
                    }
                }

                // schedule the patient
                value.scanWeek = week[i];
                value.scanDay = day[i];
                value.slotNr = slot[i];
                value.appTime = weekSchedule[day[i]][slot[i]].appTime;

                // update moving average elective appointment waiting time
                if (value.patientType == 1) {
                    if (previousWeek < week[i]) {
                        movingAvgElectiveAppWT[previousWeek] = movingAvgElectiveAppWT[previousWeek] / numberOfElectivePerWeek;
                        numberOfElectivePerWeek = 0;
                        previousWeek = week[i];
                    }
                    wt = value.getAppWT();
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
                            if (weekSchedule[d][s].patientType == value.patientType) {
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
        patients.sort((p1, p2) -> {
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
                if (p1.patientType == 2)
                    return 1;
                if (p2.patientType == 2)
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
            if (p1.patientType == 2)
                return 1;
            if (p2.patientType == 2)
                return 0;
            if (p1.nr < p2.nr)
                return 1;
            if (p1.nr > p2.nr)
                return 0;
            return 1;
        });
    }

    public double exponential_distribution(double lambda) {
        Random r = generateRandom();
        float j1 = (float) r.nextInt(1000 + 1) / 1000;
        if (j1 == 0) {
            j1 += 0.0001;
        }
        float j2 = (float) (-log(j1) / lambda);

        return j2;
    }

    public int normal_distribution(double mean, double stdev) {
        // TO MODEL BASED ON CUMULATIVE DENSITY FUNCTION OF NORMAL DISTRIBUTION BASED ON BOOK OF SHELDON ROSS, Simulation, The polar method, p80.

        float v1, v2, t;
        int x;
        Random r = generateRandom();
        do {
            //System.out.println("this is the random from norm: " + r.nextFloat());
            v1 = (float) r.nextInt(1000 + 1) * 2;
            //System.out.println("this is random number exp: " + v1);
            v1 /= 1000;
            v1 -= 1;
            v2 = (float) r.nextInt(1000 + 1) * 2;
            v2 /= 1000;
            v2 -= 1;
            t = v1 * v1 + v2 * v2;
        }
        while (t >= 1 || t == 0);
        float multiplier = (float) sqrt(-2 * log(t) / t);
        x = (int) (v1 * multiplier * stdev + mean);
        return x;
    }

    public boolean bernouilli_distribution(double prob) {
        Random r = generateRandom();
        //System.out.println("this is the random from bern: " + seed);
        float j1 = (float) r.nextInt(1000 + 1) / 1000;
        //System.out.println("this is random number exp: " + j1);
        if (j1 < prob)
            return false;
        else
            return true;
    }
}

