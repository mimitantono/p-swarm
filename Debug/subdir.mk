################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
O_SRCS += \
../cluster.o \
../cpu_info.o \
../db.o \
../main.o \
../matrix.o \
../nw.o \
../property.o \
../qgram.o \
../scan.o \
../search16.o \
../search8.o \
../util.o 

CC_SRCS += \
../cluster.cc \
../cpu_info.cc \
../db.cc \
../main.cc \
../matrix.cc \
../nw.cc \
../property.cc \
../qgram.cc \
../scan.cc \
../search16.cc \
../search8.cc \
../util.cc 

OBJS += \
./cluster.o \
./cpu_info.o \
./db.o \
./main.o \
./matrix.o \
./nw.o \
./property.o \
./qgram.o \
./scan.o \
./search16.o \
./search8.o \
./util.o 

CC_DEPS += \
./cluster.d \
./cpu_info.d \
./db.d \
./main.d \
./matrix.d \
./nw.d \
./property.d \
./qgram.d \
./scan.d \
./search16.d \
./search8.d \
./util.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -I"/Users/mimitantono/Documents/workspace/pcluster/cityhash" -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


