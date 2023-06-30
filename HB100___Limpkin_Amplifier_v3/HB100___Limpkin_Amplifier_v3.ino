#define cbi(sfr, bit) (_SFR_BYTE(sfr) &= ~_BV(bit))
#define sbi(sfr, bit) (_SFR_BYTE(sfr) |= _BV(bit))

const byte          hb100_fout_pin      =   2;     // interrupt capable pin (2 or 3 for Nano)
const byte          hb100_vout_pin      =  A0;     // analog input pin
const unsigned long analysis_every      = 500;     // ms
const float         period_to_frequency =  70.235; // Hz/m/s
const byte          level_minimum       = 120;     // minimum signal level for analysis
const float         frequency_lower     =   4.0;   // Hz
const float         frequency_upper     =  70.0;   // Hz
const byte          average_pulses      =   4;     // min 1

volatile unsigned long pulse_times[256];  // pulses start times, 256 element ringbuffer, use with counter of type byte (!)
volatile byte          pulse_levels[256]; // pulses signal levels, 256 element ringbuffer, use with counter of type byte (!)
volatile byte          pulse_last;        // last pulse written in pulse_times ring buffer

void setup()
{
   unsigned int i;

   for (i = 0; i < 256; i++) // initialize pulse_times ring buffer
   {
       pulse_times[i] = 0; 
   }
   
   pulse_last = 255; // initialize pulse_last

   // Increase ADC sample rate from 10kS/sec (100us/S) to 1MS/s (1us/S)
   sbi(ADCSRA, ADPS2);
   cbi(ADCSRA, ADPS1);
   cbi(ADCSRA, ADPS0);
   
   Serial.begin(250000); 
   pinMode(hb100_fout_pin, INPUT);
   pinMode(hb100_vout_pin, INPUT);
   attachInterrupt(digitalPinToInterrupt(hb100_fout_pin), hb100InterruptFoutPin, RISING);
}

void loop()
{
   static unsigned long analysis_next = millis();

   byte          pulse_current;
   byte          pulse_previouse;       
   unsigned long time_current;
   unsigned long time_previouse;
   byte          level_current;
   byte          level_previouse;
   long          period;
   float         frequency;
   float         velocity;
   
   if (millis() >= analysis_next)
   {
      pulse_current   = pulse_last;  // copy volatile variable
      pulse_previouse = pulse_current;
      pulse_previouse = pulse_previouse - average_pulses;

      time_current    = pulse_times[pulse_current];   
      time_previouse  = pulse_times[pulse_previouse];
      level_current   = pulse_levels[pulse_current];
      level_previouse = pulse_levels[pulse_current];
      
      period         = (time_current - time_previouse) / average_pulses;

      if (time_current / 1000 > millis() - analysis_every && 
          period              > 0                         &&
          level_current       > level_minimum             &&
          level_previouse     > level_minimum               )
      {
         frequency = 1000000.0 / period;

         //if (frequency > frequency_upper || frequency < frequency_lower)
         //{
         //   frequency = 0.0;
         //}
         
         velocity = frequency / period_to_frequency;
      }
      else
      {
         frequency = 0.0;
         velocity  = 0.0;
      }
      
      /*Serial.print(frequency, 2);*/
      /*Serial.print("\t");*/
      Serial.print(velocity * 100, 3);
      Serial.print("\t");
      Serial.print(level_current / 2, DEC);
      Serial.println("");
         
      analysis_next = millis() + analysis_every;
   }
}

void hb100InterruptFoutPin()
{
   pulse_last++;
   pulse_times[pulse_last] = micros();
   pulse_levels[pulse_last] = analogRead(hb100_vout_pin) / 4;
}
