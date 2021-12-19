// teensy 3.6 sketch to implement closed loop OMR

#include <math.h>
#include <TimerOne.h>
#include <CircularBuffer.h>

#define CIRCULAR_BUFFER_INT_SAFE

/*Frequency at which the data is acquired in samples/second*/
const long sampleFreq = 10000;

/*Frequency at which the data is written*/
//const int writeFreq = 1000;

/*buffer size in ms*/
const int bufferSz = 10;

/*buffer size in samples*/
const int bufferSamples = sampleFreq * bufferSz / 1000;

/*ring buffer the VNR data is read to*/
CircularBuffer<volatile float, bufferSamples> readBuffer;

/*EMA filter alpha*/
const float alpha = 0.005;

// Common global variables

elapsedMillis triggerTime;           // msec, time since last EPSC was triggered

int currAIval = 0; //analog in value
int currAOval = 0; //analog out value of stdev processed burst signal
int currAOval2 = 0; // analog out value of integrated speed signal
float currSum = 0; //current sum of values in circular buffer
float currMean = 0; // current mean of values in circular buffer
float currStd = 0; // current Stdev of values in circular buffer
float vFiltLowLast = 0; // last sample of low pass filtered data
float vFiltLowCurr = 0; // current sample of low pass filtered data
float vFiltHigh = 0; //current sample of high pass filtered data
float currAcc = 0; // current acceleration, area under the curve of the processed signal
float currVelocity = 0; //current speed, integration of the acceleration
float lastVelocity = 0; //speed at iteration n-1
int tau = 1000;         // time constant of speed decrease;

bool bufferReady = true; // marks if the circular buffer is not used by another process
bool newValAvail = false; // marks if a new value is available in the circular buffer
bool period = false; // for debugging

/*parameters sent by the host computer over the USB port*/
const int nPars = 1;              // number of adjustable parameters
float thresh = 100;               // threshold for burst detection on processed signal


// Hardware connections
const int analogInPin = 0;        // ADC pin used to read VNR signal
const int analogOutPin1 = A21;     // DAC pin used to output processed burst signal
const int analogOutPin2 = A22;    // DAC pin used to output calculated speed signal
const int TriggerPin = 2;     // pin used to trigger (unused)
const int burstTriggerPin = 3; // pin used to output when burst occurs
const int led = 13;

//called by interrupt
void readAItoBuffer() {
  //read AI to circular buffer

  //  period = !period;

  // check if buffer is not currently read by main loop to avoid accessing the memory at the same time
  if (bufferReady) {
    currAIval = analogRead(analogInPin);

    //EMA low pass filter using previous value

    vFiltLowCurr = (alpha * currAIval + (1 - alpha) * vFiltLowLast);
    vFiltLowLast = vFiltLowCurr;
    vFiltHigh = currAIval - vFiltLowCurr;

    currSum = currSum - readBuffer.first() + vFiltHigh;

    // push high pass filtered signal to buffer
    readBuffer.push(vFiltHigh);

    currMean = currSum / readBuffer.size();

    // signal new samples in the buffer
    newValAvail = true;
  }

}



void setup() {
  Serial.begin(115200);
  analogWriteResolution(16);
  analogReadResolution(16);
  while (Serial.available() > 0) {              // make sure serial buffer is clear
    char foo = Serial.read();
  }
  pinMode(TriggerPin, INPUT_PULLDOWN);               // define the pin that receives triggers
  pinMode(led, OUTPUT);                   // LED to indicate bursts happening
  pinMode(burstTriggerPin, OUTPUT);       // output pin to trigger equipment when bursts are happening

  /*period at 1kHz sampling is 100 us, should use 'sampleFreq', but didn't work
     hardcoded for now
  */
  Timer1.initialize(100);
  // Timer1.initialize((1 / sampleFreq) * 1000000);
  Timer1.attachInterrupt(readAItoBuffer);
  //  Timer1.start();

  //push one value to buffer before getting started
  readBuffer.push(0);
}




void loop() {

  //check Serial Communication & Parameter Setting

  union {
    byte asByte[4];
    float asFloat;
  } data1;
  static float dataTemp[nPars] = {0.0};
  static byte parNo = 0;
  if (Serial.available() > 3) {
    for (int x = 0; x < 4; x++) data1.asByte[x] = Serial.read();
    dataTemp[parNo] = data1.asFloat;
    parNo++;
    if (parNo == nPars) {
      thresh = dataTemp[0];
      parNo = 0;
    }
  }

  //check if the interrup pushed a new sample to the buffer
  if (newValAvail) {
    //period = !period;

    //block the buffer from being read by the interrupt
    bufferReady = false;
    float sumDiff = 0;
    // loop over buffer to calculated the sum of the absolute differece from the mean
    int i = 0;
    while (i < readBuffer.size()) {
      sumDiff = sumDiff + sq((readBuffer[i] - currMean));
      i++;
    }
    //signal that the latest sample has bee processed
    newValAvail = false;

    // free the buffer
    bufferReady = true;

    // calculate Stdev
    currStd = sqrt(sumDiff / readBuffer.size());

    currAOval = constrain((int)currStd * 100, 0, 65536);
    // output trigger when fish is swimming and blink LED
    if (currStd > thresh) {
      digitalWrite(led, HIGH);
      digitalWrite(burstTriggerPin, HIGH);
      // integrate values of stdev processes burst signal
      currAcc = currAcc + currStd / 10000;
    }
    else {
      currAcc = 0;
      digitalWrite(led, LOW);
    }
    // calculate speed signal
    currVelocity = currAcc + lastVelocity - lastVelocity / tau;
    lastVelocity = currVelocity;
    currAOval2 = constrain((int)currVelocity * 10, 0, 65536);

    //write swim signal and speed signal to AO
    analogWrite(analogOutPin1, currAOval);
    analogWrite(analogOutPin2, currAOval2);

  }


  /*for debugging: print AO values to serial*/
  if (triggerTime > 10) {
    //  Serial.print(currAIval);
    //  Serial.print("\t");
    Serial.print(currStd);
    Serial.print("\t");
    Serial.print(thresh);
    Serial.println();
    //period = !period;
    triggerTime = 0;
  }


  //
  //if (period) {
  //  analogWrite(analogOutPin1, 0);
  //}
  //else {
  //  analogWrite(analogOutPin1, 4095);
  //}

}

