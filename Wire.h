/*
  TwoWire.h - TWI/I2C library for Arduino & Wiring
  Copyright (c) 2006 Nicholas Zambetti.  All right reserved.
  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.
  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.
  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
  Modified 2012 by Todd Krein (todd@krein.org) to implement repeated starts
*/

#ifndef TwoWire_h
#define TwoWire_h

#include <inttypes.h>
#include "Stream.h"


#ifdef USE_USI_I2C
#define USI_SEND         0              // indicates sending to TWI
#define USI_RCVE         1              // indicates receiving from TWI
#define USI_BUF_SIZE    19              // bytes in message buffer

//class USI_TWI : public Stream
class USI_TWI
{
  private:
	static uint8_t USI_Buf[];           // holds I2C send and receive data
	static uint8_t USI_BufIdx;          // current number of bytes in the send buff
	static uint8_t USI_LastRead;        // number of bytes read so far
	static uint8_t USI_BytesAvail;      // number of bytes requested but not read
	
  public:
    USI_TWI();
    void    begin();
    void    beginTransmission(uint8_t);
    size_t  write(uint8_t);
    inline size_t write(uint8_t* d, uint8_t n) { uint16_t i; for (i = 0; i < n; i++) write(d[i]); return (size_t)n; }
    inline size_t write(unsigned long n) { return write((uint8_t)n); }
    inline size_t write(long n) { return write((uint8_t)n); }
    inline size_t write(unsigned int n) { return write((uint8_t)n); }
    inline size_t write(int n) { return write((uint8_t)n); }
    void send(uint8_t b)               { write(b); }
    void send(uint8_t *d, uint8_t n)   { write(d, n); }
    void send(int n)                   { write((uint8_t)n); }
    uint8_t endTransmission();
    uint8_t endTransmission(uint8_t);
    uint8_t requestFrom(uint8_t, uint8_t);
    int     read(); 
    int     available(); 
    int     peek(void);
    void    flush(void);
    uint8_t receive(void) {
        int c = read();
        if (c < 0) return 0;
        return c;
    }
};

extern USI_TWI TinyWireM;

#endif

#ifdef USE_SOFT_I2C


class SoftI2CMaster
{

private:
  // per object data
  uint8_t _sclPin;
  uint8_t _sdaPin;
  uint8_t _sclBitMask;
  uint8_t _sdaBitMask;
  volatile uint8_t *_sclPortReg;
  volatile uint8_t *_sdaPortReg;
  volatile uint8_t *_sclDirReg;
  volatile uint8_t *_sdaDirReg;

  uint8_t usePullups;
  
  // private methods

  void i2c_writebit( uint8_t c );
  uint8_t i2c_readbit(void);
  void i2c_init(void);
  void i2c_start(void);
  void i2c_repstart(void);
  void i2c_stop(void);
  uint8_t i2c_write( uint8_t c );
  uint8_t i2c_read( uint8_t ack );
  
public:
  // public methods
  SoftI2CMaster();
  SoftI2CMaster(uint8_t sclPin, uint8_t sdaPin);
  SoftI2CMaster(uint8_t sclPin, uint8_t sdaPin, uint8_t usePullups);

  void setPins(uint8_t sclPin, uint8_t sdaPin, uint8_t usePullups);

  uint8_t beginTransmission(uint8_t address);
  uint8_t beginTransmission(int address);
  uint8_t endTransmission(void);
  uint8_t write(uint8_t);
  void write(uint8_t*, uint8_t);
  void write(int);
  void write(char*);
  void begin(void) {return;};
  uint8_t requestFrom(int address);
  uint8_t requestFrom(uint8_t address);
  uint8_t read( uint8_t ack );
  uint8_t read();
  uint8_t readLast();

};

extern SoftI2CMaster Wire;

#endif

#ifdef USE_HARD_I2C

#define BUFFER_LENGTH 32

// WIRE_HAS_END means Wire has end()
#define WIRE_HAS_END 1

class TwoWire : public Stream
{
  private:
    static uint8_t rxBuffer[];
    static uint8_t rxBufferIndex;
    static uint8_t rxBufferLength;

    static uint8_t txAddress;
    static uint8_t txBuffer[];
    static uint8_t txBufferIndex;
    static uint8_t txBufferLength;

    static uint8_t transmitting;
    static void (*user_onRequest)(void);
    static void (*user_onReceive)(int);
    static void onRequestService(void);
    static void onReceiveService(uint8_t*, int);
  public:
    TwoWire();
    void begin();
    void begin(uint8_t);
    void begin(int);
    void end();
    void setClock(uint32_t);
    void beginTransmission(uint8_t);
    void beginTransmission(int);
    uint8_t endTransmission(void);
    uint8_t endTransmission(uint8_t);
    uint8_t requestFrom(uint8_t, uint8_t);
    uint8_t requestFrom(uint8_t, uint8_t, uint8_t);
	uint8_t requestFrom(uint8_t, uint8_t, uint32_t, uint8_t, uint8_t);
    uint8_t requestFrom(int, int);
    uint8_t requestFrom(int, int, int);
    virtual size_t write(uint8_t);
    virtual size_t write(const uint8_t *, size_t);
    virtual int available(void);
    virtual int read(void);
    virtual int peek(void);
    virtual void flush(void);
    void onReceive( void (*)(int) );
    void onRequest( void (*)(void) );

    inline size_t write(unsigned long n) { return write((uint8_t)n); }
    inline size_t write(long n) { return write((uint8_t)n); }
    inline size_t write(unsigned int n) { return write((uint8_t)n); }
    inline size_t write(int n) { return write((uint8_t)n); }
    using Print::write;
};

extern TwoWire Wire;

#endif

#endif
