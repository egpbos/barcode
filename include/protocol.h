/*
 * Barcode
 * Copyright E.G.P. Bos and F.S. Kitaura
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */


#pragma once
void INIT_PROTOCOL_CONV(struct DATA *data);

void UPDATE_PROTOCOL_CONV(ULONG it, real_prec res, struct DATA *data);

void PROTOCOL_RESTART(ULONG it, struct DATA *data);

void INIT_PROTOCOL_SPEC(struct DATA *data);

void UPDATE_PROTOCOL_SPEC(struct DATA *data, const std::string & specname);







