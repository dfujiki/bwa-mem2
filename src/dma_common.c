/*
 * Amazon FPGA Hardware Development Kit
 *
 * Copyright 2018 Amazon.com, Inc. or its affiliates. All Rights Reserved.
 *
 * Licensed under the Amazon Software License (the "License"). You may not use
 * this file except in compliance with the License. A copy of the License is
 * located at
 *
 *    http://aws.amazon.com/asl/
 *
 * or in the "license" file accompanying this file. This file is distributed on
 * an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, express or
 * implied. See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <fcntl.h>
#include <stdio.h>
#include <errno.h>
#include <stdint.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "dma_common.h"
//#include "user_defines.h"

#define AMZ_PCI_VENDOR_ID   0x1D0F
#define PCI_DEVICE_ID       0xF000

int fill_buffer_urandom(uint8_t *buf, size_t size)
{
    int fd, rc;
    off_t i = 0;

    fd = open("/dev/urandom", O_RDONLY);
    if (fd < 0) {
        return errno;
    }

    for (i = 0; i < size; ) {
        rc = read(fd, buf + i, min(4096, size - i));
        if (rc < 0) {
            close(fd);
            return errno;
        }
        i += rc;
    }
    close(fd);

    return 0;
}

uint64_t buffer_compare(uint8_t *bufa, uint8_t *bufb,
    size_t buffer_size)
{
    size_t i;
    uint64_t differ = 0;
    for (i = 0; i < buffer_size; ++i) {
        if (bufa[i] != bufb[i]) {
            differ += 1;
        }
    }

    return differ;
}

int check_slot_config(int slot_id)
{
    int rc;
    struct fpga_mgmt_image_info info = {0};

    /* get local image description, contains status, vendor id, and device id */
    rc = fpga_mgmt_describe_local_image(slot_id, &info, 0);
    if(rc != 0){
        printf("Unable to get local image information. Are you running as root?\n"); 
    }
    fail_on(rc, out, "Unable to get local image information. Are you running "
        "as root?");

    /* check to see if the slot is ready */
    if (info.status != FPGA_STATUS_LOADED) {
        rc = 1;
        printf("Slot %d is not ready",slot_id);
        fail_on(rc, out, "Slot %d is not ready", slot_id);
    }

    /* confirm that the AFI that we expect is in fact loaded */
    if (info.spec.map[FPGA_APP_PF].vendor_id != AMZ_PCI_VENDOR_ID ||
        info.spec.map[FPGA_APP_PF].device_id != PCI_DEVICE_ID)
    {
        rc = 1;
        char sdk_path_buf[512];
        char *sdk_env_var;
        sdk_env_var = getenv("SDK_DIR");
        snprintf(sdk_path_buf, sizeof(sdk_path_buf), "%s",
            (sdk_env_var != NULL) ? sdk_env_var : "<aws-fpga>");
        log_error(
            "...\n"
            "  The slot appears loaded, but the pci vendor or device ID doesn't match the\n"
            "  expected values. You may need to rescan the fpga with \n"
            "    fpga-describe-local-image -S %i -R\n"
            "  Note that rescanning can change which device file in /dev/ a FPGA will map to.\n",
            slot_id);
        log_error(
            "...\n"
            "  To remove and re-add your xdma driver and reset the device file mappings, run\n"
            "    sudo rmmod xdma && sudo insmod \"%s/sdk/linux_kernel_drivers/xdma/xdma.ko\"\n",
            sdk_path_buf);
        fail_on(rc, out, "The PCI vendor id and device of the loaded image are "
                         "not the expected values.");
    }

out:
    return rc;
}

int initialize_read_queue(int slot_id, int channel) {
    int read_fd = -1;

    read_fd = fpga_dma_open_queue(FPGA_DMA_XDMA, slot_id,
				/*channel*/ channel, /*is_read*/ true);

    if(read_fd < 0){
        printf("Unable to open read dma queue\n");
        return -1;
    }

    return read_fd;
}

int initialize_write_queue(int slot_id, int channel) {
    int write_fd = -1;

    write_fd = fpga_dma_open_queue(FPGA_DMA_XDMA, slot_id,
            /*channel*/ channel, /*is_read*/ false);

    if(write_fd < 0){
        printf("Unable to open write dma queue\n");
        return -1;
    }
    
    return write_fd;
}


void close_read_queue(int read_fd){
    close(read_fd);
    return;
}

void close_write_queue(int write_fd) {
    close(write_fd);
    return;
}

pci_bar_handle_t initialize_ocl_bus(int slot_id){
    int ret = 0;
    pci_bar_handle_t pci_bar_handle = PCI_BAR_HANDLE_INIT; 

    ret = fpga_pci_attach(slot_id,FPGA_APP_PF,APP_PF_BAR0,0,&pci_bar_handle);

    return pci_bar_handle;
}

void close_ocl_bus(pci_bar_handle_t pci_bar_handle){
   
    int ret = 0;

    ret = fpga_pci_detach(pci_bar_handle);

    return;
}



int write_to_fpga(int write_fd, uint8_t * buffer, int len, uint64_t addr){
   
    int ret = 0;
    ret = fpga_dma_burst_write(write_fd, buffer, len, addr);

    if(ret < 0) {
        fprintf(stderr, "Burst write failed on 0x%lx \n",addr);
    }
    return ret;
}

// Allocated inside the function, should be freed outside

uint8_t * read_from_fpga(int read_fd, int len, uint64_t addr){
   
    int ret = 0;
    uint8_t * buffer = (uint8_t * ) malloc (len * sizeof(uint8_t));

    ret = fpga_dma_burst_read(read_fd, buffer, len, addr);

    if(ret < 0){
        fprintf(stderr, "Burst read failed on 0x%lx \n",addr);
        return NULL;
    }

    return buffer;
}



void print_byte_hex(uint8_t data){
    uint8_t mask = 0x0f;
    uint8_t result = data >> 4;
    printf("%x",result);
    result = data & mask;
    printf("%x",result);
}

void print_byte_bin(uint8_t data){
    uint8_t mask = 0x80;
    uint8_t result = data & mask;
    int i = 0;
    for(i=0;i<8;i++){
        result = data & mask;
        if(result != 0){
            printf("1");
        }
        else{
            printf("0");
        }
        mask = mask >> 1;
    }
}

void print_data(uint8_t * data_buffer, int data_buffer_size, int format, int partition){
    int i = 0;
    for(i=data_buffer_size - 1;i>=0;i--){
        
        if(partition != 0){
            if((data_buffer_size - 1 - i) != 0){
                if(((data_buffer_size - 1 - i) % partition) == 0){
                    printf("|");
                }
            }
        }

        if(format == 1){
            print_byte_bin(data_buffer[i]);
        }
        else{
            print_byte_hex(data_buffer[i]);
        }
    }
}

void print_data_cacheblock_aligned(uint8_t * data_buffer, int data_buffer_size, int format, int partition){
    int num_cache_blocks = data_buffer_size / CACHE_BLOCK_SIZE;
    if(data_buffer_size % CACHE_BLOCK_SIZE != 0) {
        num_cache_blocks += 1;
    }
    
    int i = 0;
    int j = 0;
    int k = 0;

    int new_partition = 0;

    for(i = 0;i<num_cache_blocks;i++){
        if(partition != 0){
            new_partition = CACHE_BLOCK_SIZE % partition;
            if(new_partition == 0){
                new_partition = partition;
            }
        }
        for(j = CACHE_BLOCK_SIZE-1,k = 0;j>=0;j--, k++){
            int index = (i * CACHE_BLOCK_SIZE) + j;
            if(partition != 0 && k != 0){
                if(((CACHE_BLOCK_SIZE - k) != 0) && (((CACHE_BLOCK_SIZE - k) % new_partition) == 0)){
                    printf("|"); 
                }
                new_partition = partition;
            }
            if(format == 1){
                // Print data in binary
                if(index < data_buffer_size){
                    print_byte_bin(data_buffer[index]);
                }
                else {
                    print_byte_bin((uint8_t)0);
                }
            }
            else {
                // default is hex
                if(index < data_buffer_size){
                    print_byte_hex(data_buffer[index]);
                }
                else {
                    print_byte_hex((uint8_t)0);
                }
            }


        }
        printf("\n");
    }
}
