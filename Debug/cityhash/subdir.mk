################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
O_SRCS += \
../cityhash/city.o 

CC_SRCS += \
../cityhash/city.cc 

OBJS += \
./cityhash/city.o 

CC_DEPS += \
./cityhash/city.d 


# Each subdirectory must supply rules for building sources it contributes
cityhash/%.o: ../cityhash/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -I"/Users/mimitantono/Documents/workspace/pcluster/cityhash" -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


