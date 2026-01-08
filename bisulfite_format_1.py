#!/usr/bin/env python3
"""
Скрипт для бисульфитной конверсии - ИСПРАВЛЕННАЯ ВЕРСИЯ

Биологическая логика:
=====================
Бисульфитная конверсия: C → T, но C в CpG (CG) сохраняется как N (метилированный)

Двуцепочечная ДНК (пример):
5'-ACGTCGACGT-3'  (плюс-цепь, читаем слева направо)
3'-TGCAGCTGCA-5'  (минус-цепь, комплементарна, антипараллельна)

CpG сайты:
- На плюс-цепи: позиции где CG (читая 5'→3')
- На минус-цепи: те же позиции, но комплементарно — это GC (читая 3'→5')
  Однако если прочитать минус-цепь в её естественном 5'→3' направлении,
  то это тоже будет CG!

ВАЖНО: Бисульфитная конверсия происходит на ДЕНАТУРИРОВАННОЙ (одноцепочечной) ДНК,
поэтому каждая цепь конвертируется независимо, читая в направлении 5'→3'.

Выходные файлы:
===============
Файл 1: Исходные цепи (без конверсии)
Файл 2: Обе цепи после бисульфитной конверсии  
Файл 3: Плюс-цепь конвертированная + её комплемент (для ПЦР)
Файл 4: Минус-цепь конвертированная + её комплемент (для ПЦР)

Использование: python bisulfite_format_corrected.py input.fasta
"""

import sys
import os


def bisulfite_convert(seq_str):
    """
    Бисульфитная конверсия последовательности в направлении 5'→3'.
    
    Правила:
    - CG → NG (C в CpG заменяется на N, обозначая потенциально метилированный C)
    - C (не в CpG) → T (неметилированный C конвертируется в T)
    - Все остальные основания остаются без изменений
    
    ВАЖНО: Последовательность должна быть в направлении 5'→3'!
    """
    seq = seq_str.upper()
    result = []
    i = 0
    while i < len(seq):
        if i < len(seq) - 1 and seq[i:i+2] == 'CG':
            # CpG динуклеотид: C может быть метилирован → N
            result.append('N')
            result.append('G')
            i += 2
        elif seq[i] == 'C':
            # C вне CpG контекста → T
            result.append('T')
            i += 1
        else:
            result.append(seq[i])
            i += 1
    return ''.join(result)


def complement_sequence(seq_str):
    """
    Возвращает комплементарную последовательность (без реверса).
    A↔T, G↔C
    """
    complement_dict = {
        'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
        'N': 'N', 'R': 'Y', 'Y': 'R', 'S': 'S',
        'W': 'W', 'K': 'M', 'M': 'K', 'B': 'V',
        'D': 'H', 'H': 'D', 'V': 'B'
    }
    return ''.join([complement_dict.get(base.upper(), 'N') for base in seq_str])


def reverse_complement(seq_str):
    """
    Возвращает обратно-комплементарную последовательность.
    Это то, что нужно для получения минус-цепи в направлении 5'→3'.
    """
    return complement_sequence(seq_str)[::-1]


def read_sequences(filename):
    """Чтение последовательностей из FASTA или текстового файла."""
    sequences = []
    if filename.endswith(('.fasta', '.fa', '.fna')):
        current_header = None
        current_seq = []
        with open(filename, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_header:
                        sequences.append((current_header, ''.join(current_seq)))
                    current_header = line[1:]
                    current_seq = []
                else:
                    current_seq.append(line)
            if current_header:
                sequences.append((current_header, ''.join(current_seq)))
    else:
        with open(filename, 'r') as f:
            seq = ''.join(line.strip() for line in f if line.strip())
        sequences.append(('input_sequence', seq))
    return sequences


def create_output_files(base_name, sequences):
    """
    Создание 4 выходных файлов для каждой последовательности.
    """
    for header, original_seq in sequences:
        name = header if header != 'input_sequence' else 'sequence'
        original_seq = original_seq.upper()
        
        # =====================================================
        # ОПРЕДЕЛЕНИЕ ЦЕПЕЙ
        # =====================================================
        
        # Плюс-цепь: исходная последовательность (5'→3')
        plus_5to3 = original_seq
        
        # Минус-цепь в направлении 5'→3' = reverse complement плюс-цепи
        minus_5to3 = reverse_complement(original_seq)
        
        # Минус-цепь в направлении 3'→5' = просто комплемент (для отображения)
        minus_3to5 = complement_sequence(original_seq)
        
        # =====================================================
        # БИСУЛЬФИТНАЯ КОНВЕРСИЯ
        # Конверсия всегда применяется к цепи в направлении 5'→3'!
        # =====================================================
        
        # Конверсия плюс-цепи (уже в 5'→3')
        plus_converted_5to3 = bisulfite_convert(plus_5to3)
        
        # Конверсия минус-цепи (сначала получаем её в 5'→3', потом конвертируем)
        minus_converted_5to3 = bisulfite_convert(minus_5to3)
        
        # Минус конвертированная в направлении 3'→5' (для отображения рядом с плюс-цепью)
        minus_converted_3to5 = minus_converted_5to3[::-1]
        
        # =====================================================
        # ФАЙЛ 1: Исходные цепи (без конверсии)
        # =====================================================
        with open(f"{base_name}_1_{name}.txt", 'w') as f:
            f.write(f"# Файл 1: Исходная ДНК (без конверсии)\n")
            f.write(f"# Плюс-цепь (5'→3') и комплементарная минус-цепь (3'→5')\n\n")
            f.write(f"5'-{plus_5to3}-3'\n")
            f.write(f"   {'|' * len(plus_5to3)}\n")
            f.write(f"3'-{minus_3to5}-5'\n")
        
        # =====================================================
        # ФАЙЛ 2: Обе цепи после бисульфитной конверсии
        # =====================================================
        with open(f"{base_name}_2_{name}.txt", 'w') as f:
            f.write(f"# Файл 2: Обе цепи после бисульфитной конверсии\n")
            f.write(f"# Плюс конвертированная (5'→3') и минус конвертированная (3'→5')\n")
            f.write(f"# ВНИМАНИЕ: цепи больше НЕ комплементарны друг другу!\n\n")
            f.write(f"5'-{plus_converted_5to3}-3'  (плюс, конвертированная)\n")
            f.write(f"\n")
            f.write(f"3'-{minus_converted_3to5}-5'  (минус, конвертированная)\n")
        
        # =====================================================
        # ФАЙЛ 3: Плюс-цепь конвертированная + её комплемент (для ПЦР)
        # После ПЦР получается дуплекс из конвертированной цепи
        # =====================================================
        plus_converted_complement_3to5 = complement_sequence(plus_converted_5to3)
        
        with open(f"{base_name}_3_{name}.txt", 'w') as f:
            f.write(f"# Файл 3: Плюс-цепь конвертированная + комплемент (для ПЦР)\n")
            f.write(f"# Это ПЦР-продукт с плюс-цепи после бисульфитной конверсии\n\n")
            f.write(f"5'-{plus_converted_5to3}-3'\n")
            f.write(f"   {'|' * len(plus_converted_5to3)}\n")
            f.write(f"3'-{plus_converted_complement_3to5}-5'\n")
        
        # =====================================================
        # ФАЙЛ 4: Минус-цепь конвертированная + её комплемент (для ПЦР)
        # Минус-цепь показана в направлении 5'→3'
        # =====================================================
        minus_converted_complement_3to5 = complement_sequence(minus_converted_5to3)
        
        with open(f"{base_name}_4_{name}.txt", 'w') as f:
            f.write(f"# Файл 4: Минус-цепь конвертированная + комплемент (для ПЦР)\n")
            f.write(f"# Это ПЦР-продукт с минус-цепи после бисульфитной конверсии\n")
            f.write(f"# Минус-цепь показана в её естественном направлении 5'→3'\n\n")
            f.write(f"5'-{minus_converted_5to3}-3'\n")
            f.write(f"   {'|' * len(minus_converted_5to3)}\n")
            f.write(f"3'-{minus_converted_complement_3to5}-5'\n")


def verify_conversion(test_seq="ACGTCGACGT"):
    """
    Функция для проверки корректности конверсии.
    Запускается с флагом --test
    """
    print("=" * 60)
    print("ПРОВЕРКА БИСУЛЬФИТНОЙ КОНВЕРСИИ")
    print("=" * 60)
    
    print(f"\nИсходная последовательность: {test_seq}")
    print()
    
    # Показываем исходную двуцепочечную ДНК
    plus_5to3 = test_seq.upper()
    minus_3to5 = complement_sequence(plus_5to3)
    minus_5to3 = reverse_complement(plus_5to3)
    
    print("ИСХОДНАЯ ДНК:")
    print(f"  5'-{plus_5to3}-3'  (плюс-цепь)")
    print(f"  3'-{minus_3to5}-5'  (минус-цепь)")
    print()
    print(f"  Минус-цепь в направлении 5'→3': {minus_5to3}")
    print()
    
    # Находим CpG сайты
    print("CpG САЙТЫ:")
    cpg_positions_plus = []
    for i in range(len(plus_5to3) - 1):
        if plus_5to3[i:i+2] == 'CG':
            cpg_positions_plus.append(i)
            print(f"  Плюс-цепь: позиция {i}-{i+1} = CG")
    
    cpg_positions_minus = []
    for i in range(len(minus_5to3) - 1):
        if minus_5to3[i:i+2] == 'CG':
            cpg_positions_minus.append(i)
            print(f"  Минус-цепь (5'→3'): позиция {i}-{i+1} = CG")
    print()
    
    # Бисульфитная конверсия
    plus_converted = bisulfite_convert(plus_5to3)
    minus_converted = bisulfite_convert(minus_5to3)
    
    print("БИСУЛЬФИТНАЯ КОНВЕРСИЯ:")
    print(f"  Плюс-цепь:")
    print(f"    До:    5'-{plus_5to3}-3'")
    print(f"    После: 5'-{plus_converted}-3'")
    print()
    print(f"  Минус-цепь (в направлении 5'→3'):")
    print(f"    До:    5'-{minus_5to3}-3'")
    print(f"    После: 5'-{minus_converted}-3'")
    print()
    
    # Проверка
    print("ПРОВЕРКА:")
    
    # Подсчёт C→N (в CpG) и C→T (вне CpG)
    c_to_n_plus = plus_converted.count('N')
    c_to_n_minus = minus_converted.count('N')
    
    original_c_plus = plus_5to3.count('C')
    original_c_minus = minus_5to3.count('C')
    
    converted_t_plus = plus_converted.count('T') - plus_5to3.count('T')
    converted_t_minus = minus_converted.count('T') - minus_5to3.count('T')
    
    print(f"  Плюс-цепь: {original_c_plus} C → {c_to_n_plus} N (в CpG) + {converted_t_plus} T (вне CpG)")
    print(f"  Минус-цепь: {original_c_minus} C → {c_to_n_minus} N (в CpG) + {converted_t_minus} T (вне CpG)")
    print()
    
    # Показываем выравнивание для проверки
    print("ВЫРАВНИВАНИЕ (для визуальной проверки):")
    print(f"  Плюс оригинал:       5'-{plus_5to3}-3'")
    print(f"  Плюс конверт:        5'-{plus_converted}-3'")
    print()
    print(f"  Минус оригинал:      5'-{minus_5to3}-3'")
    print(f"  Минус конверт:       5'-{minus_converted}-3'")
    print()
    
    # Детальный разбор по позициям
    print("ДЕТАЛЬНЫЙ РАЗБОР ПЛЮС-ЦЕПИ:")
    for i, (orig, conv) in enumerate(zip(plus_5to3, plus_converted)):
        if orig != conv:
            context = ""
            if i < len(plus_5to3) - 1 and plus_5to3[i:i+2] == 'CG':
                context = "(CpG контекст)"
            else:
                context = "(вне CpG)"
            print(f"    Позиция {i}: {orig} → {conv} {context}")
    print()
    
    print("ДЕТАЛЬНЫЙ РАЗБОР МИНУС-ЦЕПИ (5'→3'):")
    for i, (orig, conv) in enumerate(zip(minus_5to3, minus_converted)):
        if orig != conv:
            context = ""
            if i < len(minus_5to3) - 1 and minus_5to3[i:i+2] == 'CG':
                context = "(CpG контекст)"
            else:
                context = "(вне CpG)"
            print(f"    Позиция {i}: {orig} → {conv} {context}")
    
    print()
    print("=" * 60)


def main():
    if len(sys.argv) < 2:
        print("Использование: python bisulfite_format_corrected.py input.fasta")
        print("           или: python bisulfite_format_corrected.py --test [последовательность]")
        sys.exit(1)
    
    if sys.argv[1] == '--test':
        test_seq = sys.argv[2] if len(sys.argv) > 2 else "ACGTCGACGT"
        verify_conversion(test_seq)
        sys.exit(0)
    
    input_file = sys.argv[1]
    
    if not os.path.exists(input_file):
        print(f"Ошибка: файл '{input_file}' не найден")
        sys.exit(1)
    
    base_name = input_file.rsplit('.', 1)[0]
    sequences = read_sequences(input_file)
    
    print(f"Обработка: {input_file}")
    print(f"Найдено последовательностей: {len(sequences)}")
    
    create_output_files(base_name, sequences)
    
    print(f"Создано 4 файла для каждой последовательности")
    print("\nФайлы:")
    print("  _1_: Исходная ДНК (без конверсии)")
    print("  _2_: Обе цепи после бисульфитной конверсии")
    print("  _3_: Плюс-цепь конвертированная + комплемент (для ПЦР)")
    print("  _4_: Минус-цепь конвертированная + комплемент (для ПЦР)")


if __name__ == "__main__":
    main()