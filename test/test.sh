#!/usr/bin/env bats

setup() {
  name="samclip"
  bats_require_minimum_version 1.5.0
  dir=$(dirname "$BATS_TEST_FILENAME")
  cd "$dir"
  bin="$dir/../$name"
  exe="$bin --ref test.fna"
}

@test "Script syntax check" {
  run -0 perl -c "$bin"
}
@test "Version" {
  run -0 $bin -v
  [[ "${lines[0]}" =~ "$name" ]]
}
@test "Help" {
  run -0 $bin -h
  [[ "$output" =~ "USAGE" ]]
}
@test "No parameters" {
  run ! $bin
}
@test "Bad option" {
  run ! $bin -Y
  [[ "$output" =~ "Unknown option" ]]
  [[ ! "$output" =~ "USAGE" ]]
}
@test "Passing a folder" {
  run ! $bin $dir
  [[ "$output" =~ "ERROR" ]]
}
@test "Empty input" {
  run ! $bin /dev/null
  [[ "$output" =~ "ERROR" ]]
}

@test "Non-existent --ref" {
  run ! $bin --ref does_not_exist.fasta
  [[ "$output" =~ "ERROR" ]]
}
@test "Empty --ref index" {
  run ! $bin --ref empty.fao
  [[ "$output" =~ "ERROR" ]]
}
@test "Empty SAM file" {
  run -0 $bin --ref empty.fai /dev/null
  [[ "$output" =~ "Found 0 seq" ]]
}
@test "Empty SAM but good --ref" {
  run -0 $exe /dev/null
  [[ "$output" =~ "Found 4 seq" ]]
}
@test "Accept STDIN" {
  run -0 $exe < test.sam
  [[ "$output" =~ "SN:contig04" ]]
  [[ "$output" =~ "Header contained 5 lines" ]]
}
@test "Accept ARGV" {
  run -0 $exe test.sam
  [[ "$output" =~ "SN:contig04" ]]
  [[ "$output" =~ "Header contained 5 lines" ]]
}
@test "PAarameter --max" {
  run -0 $exe --max 0 test.sam
  [[ "$output" =~ "SN:contig04" ]]
  [[ "$output" =~ "Header contained 5 lines" ]]
}
@test "PAarameter --invert" {
  run -0 $exe --max 0 test.sam
  [[ "$output" =~ "SN:contig04" ]]
  [[ "$output" =~ "SAM records 249" ]]
}
@test "PAarameter --progress" {
  run -0 $exe --progress 50 test.sam
  [[ "$output" =~ "SN:contig03" ]]
  [[ "$output" =~ "Processed 200 records" ]]
}
