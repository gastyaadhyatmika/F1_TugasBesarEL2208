/* EL2208 Praktikum Pemecahan Masalah dengan C 2023/2024
* Modul            : Tubes - Travelling Salesmen Problem 
* Nama (NIM)       : Prajnagastya Adhyatmika (18322005)
* Asisten (NIM)    : Emmanuella Pramudita Rumanti (13220031)
* Nama File        : main.c
* Deskripsi        : penggunaan algoritma Ant Colony Optimization (ACO) untuk menyelesaikan masalah TSP (Travelling Salesman Problem)
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <time.h>

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884197159399375105820
#endif

#define MAKS_KOTA 100
#define ALPHA 1.0
#define BETA 5.0
#define RHO 0.5
#define Q 100.0
#define ITERASI_MAKSIMAL 100

struct Kota {
    char nama[50];
    double lintang;
    double bujur;
};

// Fungsi untuk menghitung jarak antara dua kota berdasarkan koordinat lintang dan bujur
double hitungJarak(struct Kota kota1, struct Kota kota2) {
    double beda_lintang = (kota2.lintang - kota1.lintang) * M_PI / 180.0;
    double beda_bujur = (kota2.bujur - kota1.bujur) * M_PI / 180.0;
    double a = pow(sin(beda_lintang / 2), 2) + cos(kota1.lintang * M_PI / 180.0) * cos(kota2.lintang * M_PI / 180.0) * pow(sin(beda_bujur / 2), 2);
    double c = 2 * atan2(sqrt(a), sqrt(1 - a));
    return 6371 * c; // Jarak dalam kilometer (dapat diubah sesuai kebutuhan)
}

// Fungsi untuk membaca data kota dari file CSV
int bacaFile(const char *nama_file, struct Kota kota[]) {
    FILE *file = fopen(nama_file, "r");
    if (file == NULL) {
        printf("Tidak dapat membuka file %s\n", nama_file);
        return -1;
    }
    int jumlah_kota = 0;
    char baris[1024];
    while (fgets(baris, sizeof(baris), file) != NULL) {
        char *nama = strtok(baris, ",");
        char *lintang_str = strtok(NULL, ",");
        char *bujur_str = strtok(NULL, ",");
        if (nama != NULL && lintang_str != NULL && bujur_str != NULL) {
            strcpy(kota[jumlah_kota].nama, nama);
            kota[jumlah_kota].lintang = atof(lintang_str);
            kota[jumlah_kota].bujur = atof(bujur_str);
            jumlah_kota++;
        }
    }
    fclose(file);
    return jumlah_kota;
}

// Fungsi untuk menghitung probabilitas pilihan berikutnya
double hitungProbabilitas(int kota_i, int kota_j, double matriks_feromon[MAKS_KOTA][MAKS_KOTA], double matriks_jarak[MAKS_KOTA][MAKS_KOTA], bool dikunjungi[MAKS_KOTA], int jumlah_kota) {
    double feromon = matriks_feromon[kota_i][kota_j];
    double jarak = matriks_jarak[kota_i][kota_j];
    double visibilitas = 1.0 / jarak; // Sederhana, bisa ditingkatkan
    if (!dikunjungi[kota_j]) {
        return pow(feromon, ALPHA) * pow(visibilitas, BETA);
    } else {
        return 0.0;
    }
}

// Fungsi untuk memilih kota berikutnya berdasarkan probabilitas
int nextCity(int kota_sekarang, double matriks_feromon[MAKS_KOTA][MAKS_KOTA], double matriks_jarak[MAKS_KOTA][MAKS_KOTA], bool dikunjungi[MAKS_KOTA], int jumlah_kota) {
    double total_probabilitas = 0.0;
    double probabilitas[MAKS_KOTA];
    for (int i = 0; i < jumlah_kota; i++) {
        if (!dikunjungi[i]) {
            probabilitas[i] = hitungProbabilitas(kota_sekarang, i, matriks_feromon, matriks_jarak, dikunjungi, jumlah_kota);
            total_probabilitas += probabilitas[i];
        } else {
            probabilitas[i] = 0.0;
        }
    }
    double nilai_acak = (double)rand() / RAND_MAX;
    double kumulatif_probabilitas = 0.0;
    for (int i = 0; i < jumlah_kota; i++) {
        if (!dikunjungi[i]) {
            kumulatif_probabilitas += probabilitas[i] / total_probabilitas;
            if (nilai_acak <= kumulatif_probabilitas) {
                return i;
            }
        }
    }
    // Jika tidak ada kota yang dipilih, kembalikan kota yang belum dikunjungi pertama kali
    for (int i = 0; i < jumlah_kota; i++) {
        if (!dikunjungi[i]) {
            return i;
        }
    }
    return -1; // Error handling jika semua kota sudah dikunjungi
}

// Fungsi untuk memperbarui jejak feromon setelah setiap semut menyelesaikan tur
void renewFeromon(double matriks_feromon[MAKS_KOTA][MAKS_KOTA], double delta_feromon[MAKS_KOTA][MAKS_KOTA], int tur[], double panjang_tur, int jumlah_kota) {
    for (int i = 0; i < jumlah_kota; i++) {
        int dari = tur[i];
        int ke = tur[(i + 1) % jumlah_kota];
        delta_feromon[dari][ke] += Q / panjang_tur;
        delta_feromon[ke][dari] += Q / panjang_tur;
    }
    for (int i = 0; i < jumlah_kota; i++) {
        for (int j = 0; j < jumlah_kota; j++) {
            matriks_feromon[i][j] = matriks_feromon[i][j] * (1.0 - RHO) + delta_feromon[i][j];
            delta_feromon[i][j] = 0.0;
        }
    }
}

// Fungsi untuk melakukan tur semut
void turSemut(double matriks_feromon[MAKS_KOTA][MAKS_KOTA], double matriks_jarak[MAKS_KOTA][MAKS_KOTA], int tur[], int jumlah_kota, int kota_awal) {
    bool dikunjungi[MAKS_KOTA];
    double delta_feromon[MAKS_KOTA][MAKS_KOTA];
    for (int i = 0; i < jumlah_kota; i++) {
        dikunjungi[i] = false;
        for (int j = 0; j < jumlah_kota; j++) {
            delta_feromon[i][j] = 0.0;
        }
    }
    int kota_sekarang = kota_awal; // Pilih kota awal berdasarkan input pengguna
    dikunjungi[kota_sekarang] = true;
    tur[0] = kota_sekarang;
    for (int i = 1; i < jumlah_kota; i++) {
        int kota_berikutnya = nextCity(kota_sekarang, matriks_feromon, matriks_jarak, dikunjungi, jumlah_kota);
        tur[i] = kota_berikutnya;
        dikunjungi[kota_berikutnya] = true;
        kota_sekarang = kota_berikutnya;
    }
}

// Fungsi untuk inisialisasi matriks jarak antara kota-kota
void inisialisasiMatriksJarak(struct Kota kota[], int jumlah_kota, double matriks_jarak[MAKS_KOTA][MAKS_KOTA]) {
    for (int i = 0; i < jumlah_kota; i++) {
        for (int j = 0; j < jumlah_kota; j++) {
            if (i == j) {
                matriks_jarak[i][j] = 0.0;
            } else {
                matriks_jarak[i][j] = hitungJarak(kota[i], kota[j]);
            }
        }
    }
}

// Fungsi utama
int main() {
    srand(time(NULL)); // Inisialisasi seed untuk angka acak
    struct Kota kota[MAKS_KOTA];
    double matriks_jarak[MAKS_KOTA][MAKS_KOTA];
    double matriks_feromon[MAKS_KOTA][MAKS_KOTA];
    int tur_terbaik[MAKS_KOTA];
    double panjang_tur_terbaik = 99999999.0; // Nilai besar untuk inisialisasi
    char nama_file[100];
    
    printf("Masukkan nama file: ");
    scanf("%s", nama_file);
    int jumlah_kota = bacaFile(nama_file, kota);
    if (jumlah_kota == -1) {
        return -1;
    }

    char nama_kota_awal[50];
    printf("Masukkan kota keberangkatan: ");
    scanf("%s", nama_kota_awal);
    int kota_awal = -1;
    for (int i = 0; i < jumlah_kota; i++) {
        if (strcmp(kota[i].nama, nama_kota_awal) == 0) {
            kota_awal = i;
            break;
        }
    }
    if (kota_awal == -1) {
        printf("Kota %s tidak ditemukan.\n", nama_kota_awal);
        return -1;
    }

    inisialisasiMatriksJarak(kota, jumlah_kota, matriks_jarak);
    // Inisialisasi matriks feromon dengan nilai awal
    for (int i = 0; i < jumlah_kota; i++) {
        for (int j = 0; j < jumlah_kota; j++) {
            matriks_feromon[i][j] = 0.1;
        }
    }

    // Mulai pengukuran waktu sebelum algoritma ACO
    clock_t waktu_mulai = clock();

    for (int iter = 0; iter < ITERASI_MAKSIMAL; iter++) {
        // Lakukan tur untuk setiap semut
        for (int semut = 0; semut < jumlah_kota; semut++) {
            int tur[MAKS_KOTA];
            turSemut(matriks_feromon, matriks_jarak, tur, jumlah_kota, kota_awal);
            double panjang_tur = 0.0;
            for (int i = 0; i < jumlah_kota; i++) {
                panjang_tur += matriks_jarak[tur[i]][tur[(i + 1) % jumlah_kota]];
            }
            if (panjang_tur < panjang_tur_terbaik) {
                panjang_tur_terbaik = panjang_tur;
                for (int i = 0; i < jumlah_kota; i++) {
                    tur_terbaik[i] = tur[i];
                }
            }
        }
        // Perbarui jejak feromon
        double delta_feromon[MAKS_KOTA][MAKS_KOTA] = {0.0};
        renewFeromon(matriks_feromon, delta_feromon, tur_terbaik, panjang_tur_terbaik, jumlah_kota);
    }

    // Akhiri pengukuran waktu setelah algoritma ACO selesai
    clock_t waktu_selesai = clock();
    double waktu_berjalan = (double)(waktu_selesai - waktu_mulai) / CLOCKS_PER_SEC;

    // Output jalur terbaik
    printf("Tur terbaik:\n");
    for (int i = 0; i < jumlah_kota; i++) {
        printf("%s -> ", kota[tur_terbaik[i]].nama);
    }
    printf("%s\n", kota[tur_terbaik[0]].nama); // Kembali ke kota awal
    printf("Jarak tur terbaik: %lf\n", panjang_tur_terbaik);
    printf("Waktu: %f detik\n", waktu_berjalan);

    return 0;
}
