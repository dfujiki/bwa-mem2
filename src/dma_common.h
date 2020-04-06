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

#pragma once

#define CACHE_BLOCK_SIZE 64

#define	MEM_16G              (1ULL << 34)

#include "fpga_pci.h"
#include "fpga_mgmt.h"
#include "fpga_dma.h"
#include "utils/lcd.h"

#ifdef __cplusplus
extern "C" {
#endif
/**
 * Fills the buffer with bytes read from urandom.
 */
int fill_buffer_urandom(uint8_t *buf, size_t size);

/**
 * This function is like memcmp, but it returns the number of bytes that differ.
 *
 * @returns number of bytes which differ, i.e. zero if buffers are the same
 */
uint64_t buffer_compare(uint8_t *bufa, uint8_t *bufb,
    size_t buffer_size);

/**
 * Checks to make sure that the slot has a recognized AFI loaded.
 */
int check_slot_config(int slot_id);

void print_data(uint8_t * data_buffer, int data_buffer_size, int format, int partition);

void print_data_cacheblock_aligned(uint8_t * data_buffer, int data_buffer_size, int format, int partition);

int initialize_read_queue(int slot_id, int channel);

int initialize_write_queue(int slot_id, int channel);

void close_read_queue(int read_fd);

void close_write_queue(int write_fd);

pci_bar_handle_t initialize_ocl_bus(int slot_id);

void close_ocl_bus(pci_bar_handle_t pci_bar_handle);

int write_to_fpga(int write_fd, uint8_t * buffer, int len, uint64_t addr);

uint8_t * read_from_fpga(int read_fd, int len, uint64_t addr);

typedef struct {
    int write_fd;
    int read_fd;
    pci_bar_handle_t pci_bar_handle;
} fpga_pci_data_t;

#ifdef __cplusplus
}
#endif