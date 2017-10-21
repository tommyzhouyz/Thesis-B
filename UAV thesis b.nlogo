extensions [array]
breed [uavs uav]
breed [targets target]
breed [hangars hangar]
breed [helis heli]

globals[
  ;Heterogeneous
  genetic-code-list
  generation-times
  global-gene-num
  globalarray
  generation-count
  GA-mode
  testarray
  homex
  homey
  MAXenergy
  MAXCount
  GBASE
  HUN
  SCM
  MAX-p
  recharge-speed
  typ1count
  typ2count
  num-genes
  gene-length
  type-gene-num
  gene-code-length
  u-avoid-u-coef
  u-avoid-t-coef
  best-pop-gene
  p-dec
]

turtles-own [
  uavs-in-range
  uavs-close-range
  targets-in-range
  targets-close-range
  force-x
  force-y
  types
]
patches-own[
p-value
]
uavs-own [
  current-target
  UU-attraction
  UU-repulsion
  UT-attraction
  UP-attraction
  p-coef
  p-inc

  UU-attraction2
  UU-repulsion2

  UU-attraction3
  UU-repulsion3

  energy
  recharging ; flag for indicating UAV recharging or not
  nearest-home
  home-dist
  uspeed

  sensor-range
]
helis-own [
  current-target
  sensor-range
]
targets-own [
  marked
  reset-counter
]
;================================
; SETUP Procedures
;================================

to setup
  clear-all
  reset-ticks
  set type-gene-num 11 ; number of genes for each type
  set testarray 10000
  if not Heterogeneous? [set type-gene-num 7]
  set global-gene-num 3
  set num-genes global-gene-num + type-gene-num * 3 ; total number of genes [2 + type-gene-num * 3]
  if not Heterogeneous? [set num-genes type-gene-num]
  set gene-length 5 ; the length in bits for each gene
  set gene-code-length gene-length * num-genes
  output-type "number of genes:"
  output-print num-genes
  ;set genetic-code-list zero-genepool ; where the initial genes are set
  set genetic-code-list random-genepool
  set generation-times n-values population-size [[]]
  ;show generation-times
  ;set globalarray n-values max-generation [[]]
  set generation-count 0
  set homex 0
  set homey 0
  set GBASE 2 ^ gene-length - 1
  set HUN 2 ^ gene-length + 0; denominator of types coef [2^gene-length]
  set SCM 2 ^ gene-length + 1 ; denominator of sensor range
  set MAX-p 255
  output-type "HUN:"
  output-print HUN
  output-type "SCM:"
  output-print SCM
  set MAXenergy 200
  set MAXCount 5 ; time for targets to reset its state after being lost
  set recharge-speed 10
  ;show genetic-code-list
  show generation-times

  ;let testlist binary-to-decimal(graycode-to-binary(item 1 genetic-code-list))

  set u-avoid-u-coef 0.6 ; the coef for sensor range to avoid uav range , must <1
  set u-avoid-t-coef 0.7 ; the coef for sensor range to avoid target range , must <1
end

;to step [gene-list]
to go
  set GA-mode 0
  ifelse generation-count <= max-generation

  [ set generation-times n-values population-size [[]]
    foreach generation-times [
      [i] -> set generation-times replace-item (position i generation-times) generation-times
             avg-time(item (position i generation-times) genetic-code-list)
      ;output-type "wut1"
      show generation-times
    ;output-type "i:"
    ;output-print i
    ]
      global-max
      ;output-type "wut"
      show generation-count
      my-update-plots
  ]
  [ stop ]
end
to go-GA
  set GA-mode 1
  ifelse generation-count <= max-generation

  [ set generation-times n-values population-size [[]]
    foreach generation-times [
      [i] -> set generation-times replace-item (position i generation-times) generation-times
             GA-fit(item (position i generation-times) genetic-code-list)
      ;output-type "wut1"
      show generation-times
    ;output-type "i:"
    ;output-print i
    ]
      global-max
      ;output-type "wut"
      show generation-count
      my-update-plots
      if testarray <= 0.05 [stop]
  ]
  [ stop ]
end

;=================================
; GA Function
;=================================
to-report GA-fit [bits]
  let zeros length filter [[i] -> i = 0] bits
  report zeros / gene-code-length
end
to-report random-genepool
  ;; local variables
  ;let gene-length 4
  ;let num-genes 29
  let chromosome-length gene-length * num-genes
  let random-code-list n-values population-size [n-values chromosome-length [one-of[0 1]]]
  ;output-type "random gene"
  ;show binary-to-decimal(graycode-to-binary(position 1 random-code-list))
  report random-code-list
end

to-report zero-genepool
  ;; local variables
  ;let gene-length 4
  ;let num-genes 29
  let chromosome-length gene-length * num-genes
  let random-code-list n-values population-size [n-values chromosome-length [0]]
  ;output-type "random gene"
  ;show binary-to-decimal(graycode-to-binary(position 1 random-code-list))
  report random-code-list
end

to-report avg-time [genetic-code] ;where the simulation runs
  let total-time 0
  let this-time 0
  ;show genetic-code
  repeat num-tests[
    setup-simulation (genetic-code)
    set this-time run-simulation
    ;show this-time
    set total-time (total-time + this-time)
  ]
  let average-time(total-time / num-tests)
  output-type "time "
  output-print average-time
  report average-time
end

to global-max

  create-next-generation
  set generation-count generation-count + 1


  output-type "time per population"
  output-print generation-times
  ;set globalarray generation-times
  let besttime min generation-times

  ifelse globalarray != []
  [if besttime <= testarray [ set testarray besttime ]]
  [set globalarray generation-times
    set testarray min generation-times]

  output-type " Best time so far: "
  output-print besttime
  output-type " Global best time so far: "
  output-print testarray
  output-type " Best genetic code "
  output-print genetic-code-list
  set best-pop-gene item 0 genetic-code-list
  output-type " Generation "
  output-print generation-count
end

to create-next-generation
  let new-generation (n-values population-size [[]])

  let tournament-pool-A-time (sublist generation-times 0 ( population-size / 2 ))
  let tournament-pool-B-time (sublist generation-times ( population-size / 2 ) population-size)
  let A-gene (sublist genetic-code-list 0 ( population-size / 2 ))
  let B-gene (sublist genetic-code-list ( population-size / 2 ) population-size)
  let parent-code-1 (item (index-of-min(tournament-pool-A-time)) A-gene)
  let parent-code-2 (item (index-of-min(tournament-pool-B-time)) B-gene)

  let new-genetic-code (crossover-recombination (parent-code-1) (parent-code-2))
  ;output-type "new-g-c:"
  ;output-print new-genetic-code
  foreach n-values (population-size / 2) [[]] [
    [i] -> let child-code-1 (item 0 new-genetic-code)
           let child-code-2 (item 1 new-genetic-code)

           set child-code-1 (random-mutation (child-code-1))
           set child-code-2 (random-mutation (child-code-2))

           ;show child-code-1
           ;show child-code-2
           let pos1 position i new-generation
           let pos2 pos1 + population-size / 2

           set new-generation replace-item pos1 new-generation child-code-1
           set new-generation replace-item pos2 new-generation child-code-2
           show new-generation
  ]

  ;let random-gene random-genepool ; add some random gene to next gen
  ;let child-code-1 (item 0 random-gene)
  ;let child-code-2 (item 1 random-gene)
  let child-code-1 parent-code-1
  let child-code-2 parent-code-2

  set new-generation replace-item 0 new-generation child-code-1
  set new-generation replace-item (0 + population-size / 2) new-generation child-code-2
  set genetic-code-list new-generation

end

to-report crossover-recombination [bits1 bits2]
  let splitsum n-values num-genes [[i] -> i * 4]
  let split-points n-of 2 splitsum

  let split-point-1 item 0 split-points
  let split-point-2 item 1 split-points
  let split-point 1 + random (length bits1 - 1)
  report list (sentence (sublist bits1 0 split-point-1)
                        (sublist bits2 split-point-1 split-point-2)
                        (sublist bits1 split-point-2 length bits1))
              (sentence (sublist bits2 0 split-point-1)
                        (sublist bits1 split-point-1 split-point-2)
                        (sublist bits2 split-point-2 length bits2))
end

to-report random-mutation [bits]
  ;let mutant bits
  set bits map [
      [b] -> ifelse-value (random-float 100 < mutation-rate)
             [1 - b]
             [b]
  ] bits
  report bits
end

;======================================
; setup Simulation functions
;======================================

to setup-simulation [genetic-code]
  clear-turtles
  reset-ticks
  let gl gene-length
  if show-background? [import-drawing "map.JPG"] ; display background image
  let type1N (binary-to-decimal(graycode-to-binary(sublist genetic-code 0 gl)))
  let type2N (binary-to-decimal(graycode-to-binary(sublist genetic-code gl (2 * gl))))
  set p-dec (binary-to-decimal(graycode-to-binary(sublist genetic-code (2 * gl) (3 * gl))))
  ;read UAV properties from genetic code
  set p-dec p-dec * MAX-pDEC / HUN
  set typ1count type1N * agent-population / HUN
  ifelse typ1count < 1 [set typ1count 1] ; bound the pop of type1 and type2
  [
    if typ1count > agent-population - 2 [set typ1count agent-population - 2]
  ]

  let leftpop agent-population - typ1count
  set typ2count type2N * agent-population / HUN
  ifelse typ2count < 1 [set typ2count 1]
  [
    if typ2count > leftpop - 1 [set typ2count leftpop - 1]
  ]
  if not Heterogeneous? [set typ1count agent-population]
  ;show genetic-code
  ;show sublist genetic-code 12 16

 ; show typ1count
  ;show typ2count
  ; create uavs
  create-uavs agent-population [
    ;set testcount testcount + 1

    let tglist type-gene-num * gene-length
    let gglist global-gene-num * gene-length
    set shape "default"
    set size 1
    set color violet
; decide type
    ifelse typ1count > 0.5
    [ set typ1count typ1count - 1
      set types 1]
    [ifelse typ2count > 0.5
      [set typ2count typ2count - 1
       set types 2]
      [set types 3]
    ]
    ifelse types = 1
    [set color yellow
      ifelse Heterogeneous? [assignprop(sublist genetic-code gglist (gglist + tglist))] ; 8 48
      [assignprop(sublist genetic-code 0 tglist)]
    ]
    [ifelse types = 2
    [set color orange
    assignprop(sublist genetic-code (gglist + tglist) (gglist + 2 * tglist))]
    [assignprop(sublist genetic-code (gglist + 2 * tglist) (gglist + 3 * tglist))]
    ]
    set energy MAXenergy

    set recharging 0
    setxy random-xcor / max-pxcor random-ycor / max-pycor
    set current-target nobody
   ; show UT-attraction
    ;show uspeed
    ;show sensor-range
  ]

  ;create targets
  create-targets target-population [
    set shape "person"
    set color blue
    set size 1
    setxy random-xcor random-ycor
    set marked 0
  ]

  create-helis floor agent-population / 2 + 1[
    set current-target nobody
    set shape "airplane"
    set color green
    set size 1
    setxy homex homey
    set sensor-range MAX-sensor-range
  ]

  create-hangars 1[
    set shape "house"
    set color red
    set size 3
    setxy homex homey
  ]
end

;======================================
; run Simulation functions
;======================================

to-report run-simulation
  ask patches [ set p-value 0 set pcolor black]
  loop[
    ask uavs[

      set force-x 0
      set force-y 0
      ifelse show-energy [set label round energy] [set label ""]
      find-targets-in-range
      find-targets-close-range
      find-uavs-in-range
      find-uavs-close-range
      set nearest-home min-one-of hangars [distance myself]
      set home-dist distance nearest-home
      if Bayesian-method?
      [
        let pinc p-inc * MAX-pINC / (2 ^ gene-length)
        ask patches in-radius sensor-range [set p-value p-value + pinc]
      ]
      ifelse recharging = 0
        [set energy energy - uspeed
          if current-target != nobody [let tsr sensor-range ask current-target [ask patches in-radius tsr [set p-value 0 set pcolor black]]]
          ifelse (energy  < home-dist)
        [
          uav-gohome

        if home-dist < 1 [set recharging 1]] ; go home if low
        ; do general stuff
        [if any? uavs-close-range    [uav-avoid-uav]           ; UAvs are repelled by other UAVs in close range
          ifelse current-target != nobody
          [uav-chase-target


          ]
          [if Bayesian-method? [uav-search-Lp]]
        ]

        ]
        [set energy energy + recharge-speed
        if energy > MAXenergy [set recharging 0]

        ]

      ifelse any? targets-in-range
      [set current-target min-one-of targets-in-range [distance myself]
        ask targets-in-range [set color red set reset-counter MAXCount]
      ]
      [ ifelse any? uavs-in-range
          [set current-target [current-target] of (min-one-of uavs-in-range [distance myself]) ]
          [set current-target nobody]
      ]

      ;if any? targets-close-range [
        ;uav-avoid-target(sensor-range)
      ;]        ; UAVs are repelled by targets in close range
      if any? uavs-in-range       [uav-cohere]              ; UAVs are attracted to other UAVs who are pursuing targets

if current-target = nobody [uav-avoid-wall]                                ; UAVs are repelled by walls in close range



    turn-towards vector-summation uav-max-turn            ; the total "force" vector determines which direction the UAV turn


    ;ask patch-ahead sensor-range [ set p-value p-value + 8 ]

    ] ; end of uavs


    ask patches with [p-value > 0] [
      if not test-bayesian? [set p-value p-value - p-dec]
      if p-value > MAX-p [set p-value MAX-p]
      if p-value < 0 [set p-value 0]
      let p-c p-value * 1.5 / 2
      set pcolor rgb p-c p-c p-c
      if not Bayesian-method? [set p-value 0 set pcolor black]
    ]
;end of patches

    ask helis[
      ;find-targets-in-range
      ;find-targets-close-range
      ask targets with [ color = red ] in-radius (2) [die]
      set targets-close-range targets with [color = red and marked = 0]
      ifelse current-target != nobody [ifelse [color] of current-target = blue [set current-target nobody][face current-target]]
      [facexy homex homey]

      if any? targets-close-range
      [ ;face min-one-of targets-close-range [distance myself]
        ;fd helis-speed
        ;ask min-one-of targets-close-range [distance myself] [set marked 1]

         ]
       ;fd helis-speed

        ;stop

    ]

    ask targets [
      set force-x 0
      set force-y 0
      if color = blue [set marked 0]
      target-find-uavs-in-range
      target-find-targets-close-range

      if any? uavs-in-range       [ target-avoid-uav ]      ; targets are repelled by UAVs within sensor range
      if any? targets-close-range [ target-avoid-target ]   ; targets are repelled by other targets in close range

      target-avoid-wall                                     ; the total "force" vector determines which direction the target turns

      turn-towards vector-summation target-max-turn
      if color = red [

        ifelse reset-counter < 0 [set color blue set marked 0]

         [set reset-counter reset-counter - 1

          if marked = 0 [
            let nearest-helis helis with [current-target = nobody]
            set nearest-helis min-one-of nearest-helis [distance myself] ; assign the closet helis to this target

            if nearest-helis != nobody [ask nearest-helis [set current-target myself] set marked 1]
          ]
         ]
      ]
      if count uavs with [ current-target = myself ] >= 1 [
        ;ask (uavs with [current-target = myself]) [set current-target nobody]
        set color red
        set reset-counter MAXCount
      ]
    ]

    ask uavs [ fd uspeed * (1 - recharging) ]
    ask targets [ fd target-speed ]
    ask helis [ fd helis-speed ]
    display
    tick

    if (count targets = 0) or (ticks >= max-run-time ) [ report ticks + (count targets) ]
  ]
end

;======================================
; uav behaviour
;======================================

to find-targets-in-range
  set targets-in-range targets in-radius sensor-range
end

to find-targets-close-range
  set targets-close-range targets in-radius (u-avoid-t-coef * sensor-range)
end

to find-uavs-in-range
  set uavs-in-range other uavs in-radius sensor-range with [current-target != nobody  ]
end

to find-uavs-close-range
  set uavs-close-range other uavs in-radius (u-avoid-u-coef * sensor-range) with [recharging = 0]
end

to uav-avoid-target [srange]
  set force-x
      force-x + UP-attraction * sum [(1 - (distance myself) / (u-avoid-t-coef * srange)) ^ 2 * sin (towards myself)] of targets-close-range
  set force-y
      force-y + UP-attraction * sum [(1 - (distance myself) / (u-avoid-t-coef * srange)) ^ 2 * cos (towards myself)] of targets-close-range
end

to uav-cohere
  if uavs-in-range with [types = 1] != nobody[
  set force-x
      force-x + UU-attraction * sum [(distance myself - 0.5 * sensor-range) / (sensor-range - 0.5 * sensor-range)
      * sin (towards myself + 180)] of uavs-in-range
  set force-y
      force-y + UU-attraction * sum [(distance myself - 0.5 * sensor-range) / (sensor-range - 0.5 * sensor-range)
      * cos (towards myself + 180)] of uavs-in-range
  ]
  if uavs-in-range with [types = 2] != nobody[
  set force-x
      force-x + UU-attraction2 * sum [(distance myself - 0.5 * sensor-range) / (sensor-range - 0.5 * sensor-range)
      * sin (towards myself + 180)] of uavs-in-range
  set force-y
      force-y + UU-attraction2 * sum [(distance myself - 0.5 * sensor-range) / (sensor-range - 0.5 * sensor-range)
      * cos (towards myself + 180)] of uavs-in-range
  ]
  if uavs-in-range with [types = 3] != nobody[
  set force-x
      force-x + UU-attraction3 * sum [(distance myself - 0.5 * sensor-range) / (sensor-range - 0.5 * sensor-range)
      * sin (towards myself + 180)] of uavs-in-range
  set force-y
      force-y + UU-attraction3 * sum [(distance myself - 0.5 * sensor-range) / (sensor-range - 0.5 * sensor-range)
      * cos (towards myself + 180)] of uavs-in-range
  ]
end
to uav-gohome
  set force-x force-x + 50 * [sin (towards myself + 180)] of nearest-home
  set force-y force-y + 50 * [cos (towards myself + 180)] of nearest-home
end
to uav-avoid-wall
  let uav-wall-separation sensor-range
  set uav-wall-separation 0.7 * uav-wall-separation
  if xcor >= (max-pxcor - uav-wall-separation)
  [set force-x force-x - (1 - (max-pxcor - uav-wall-separation) / uav-wall-separation) ^ 2]
  if xcor <= (min-pxcor + uav-wall-separation)
  [set force-x force-x + (1 - (min-pxcor + uav-wall-separation) / uav-wall-separation) ^ 2]
  if ycor >= (max-pycor - uav-wall-separation)
  [set force-y force-y - (1 - (max-pycor - uav-wall-separation) / uav-wall-separation) ^ 2]
  if ycor <= (min-pycor + uav-wall-separation)
  [set force-y force-y + (1 - (min-pycor + uav-wall-separation) / uav-wall-separation) ^ 2]
end


to uav-search-Lp ; go to areas with low p-values
  let uav-p-separation 1 * sensor-range ; could be decided by genes
  let miu 0.1
  let angle-list n-values 36 [[k] -> 5 * k]
  ;set uav-p-separation 0.8 * uav-p-separation
  let min-p MAX-p + 10000
  let min-patch nobody
  foreach angle-list ; finds the patch with the lowest p-value in the senor circle
  [[i] ->

    let temp-patch (patch-right-and-ahead i uav-p-separation)
    let t-p MAX-p + 10000
    if temp-patch != nobody
    [set t-p miu * i + [p-value] of temp-patch
    if t-p <= min-p [set min-p t-p set min-patch temp-patch]
    ]
    set temp-patch (patch-left-and-ahead i uav-p-separation)
    if temp-patch != nobody
    [set t-p miu * i + [p-value] of temp-patch
    if t-p <= min-p [set min-p t-p set min-patch temp-patch]
    ]


  ]


  let p-ahead min-patch
  ;let p-coef 0.8

  if p-ahead != nobody [
    let p-attraction (MAX-p - ([p-value] of p-ahead)) * p-coef
    set force-x force-x + p-attraction * [sin (towards myself + 180)] of p-ahead
    set force-y force-y + p-attraction * [cos (towards myself + 180)] of p-ahead
  ]

end
to uav-avoid-uav

  if uavs-close-range with [types = 1] != nobody[
  set force-x
      force-x + UU-repulsion * sum[(distance myself - u-avoid-u-coef * sensor-range) ^ 2
      * sin (towards myself) ] of uavs-close-range
  set force-y
      force-y + UU-repulsion * sum[(distance myself - u-avoid-u-coef * sensor-range) ^ 2
      * cos (towards myself) ] of uavs-close-range
  ]
  if uavs-close-range with [types = 2] != nobody[
  set force-x
      force-x + UU-repulsion2 * sum[(distance myself - u-avoid-u-coef * sensor-range) ^ 2
      * sin (towards myself) ] of uavs-close-range
  set force-y
      force-y + UU-repulsion2 * sum[(distance myself - u-avoid-u-coef * sensor-range) ^ 2
      * cos (towards myself) ] of uavs-close-range
  ]
  if uavs-close-range with [types = 3] != nobody[
  set force-x
      force-x + UU-repulsion3 * sum[(distance myself - u-avoid-u-coef * sensor-range) ^ 2
      * sin (towards myself) ] of uavs-close-range
  set force-y
      force-y + UU-repulsion3 * sum[(distance myself - u-avoid-u-coef * sensor-range) ^ 2
      * cos (towards myself) ] of uavs-close-range
  ]
end

to uav-chase-target ;; UAV procedure
  set force-x force-x + UT-attraction * [sin (towards myself + 180)] of current-target
  set force-y force-y + UT-attraction * [cos (towards myself + 180)] of current-target
end
;========================================
; UAV movements
;========================================

to turn-towards [new-heading max-turn]
  turn-at-most (subtract-headings new-heading heading) max-turn
end

to-report vector-summation
  ifelse force-x = 0 and force-y = 0
  [ report heading ]
  [ report atan force-x force-y ]
end

to turn-at-most [turn max-turn] ;; agent procedure
  ifelse abs turn > max-turn
  [ ifelse turn > 0
    [ right max-turn ]
    [ left max-turn ] ]
  [ right turn ]
end

;;=========================
;; Target Behaviours
;;=========================

to target-find-targets-close-range ;; target procedure
  set targets-close-range other targets in-radius target-target-separation
end

to target-find-uavs-in-range ;; target procedure
  set uavs-in-range uavs in-radius target-sensor-range
end

to target-avoid-uav ;; target procedure
  set force-x force-x + target-uav-repulsion * sum [(1 - (distance myself) / target-sensor-range) ^
    2 * sin (towards myself)] of uavs-in-range
  set force-y force-y + target-uav-repulsion * sum [(1 - (distance myself) / target-sensor-range) ^
    2 * cos (towards myself)] of uavs-in-range
end

to target-avoid-target ;; target procedure
  set force-x force-x + target-target-repulsion * sum [(1 - (distance myself) / target-target-separation)
    ^ 2 * sin (towards myself)] of targets-close-range
  set force-y force-y + target-target-repulsion * sum [(1 - (distance myself) / target-target-separation)
    ^ 2 * cos (towards myself)] of targets-close-range
end

to target-avoid-wall ;; target procedure
  if xcor >= (max-pxcor - target-sensor-range)
  [set force-x force-x - (1 - (max-pxcor - target-sensor-range) / target-sensor-range) ^ 2]
  if xcor <= (min-pxcor + target-target-separation)
  [set force-x force-x + (1 - (min-pxcor + target-sensor-range) / target-sensor-range) ^ 2]
  if ycor >= (max-pycor - target-target-separation)
  [set force-y force-y - (1 - (max-pycor - target-sensor-range) / target-sensor-range) ^ 2]
  if ycor <= (min-pycor + target-target-separation)
  [set force-y force-y + (1 - (min-pycor + target-sensor-range) / target-sensor-range) ^ 2]
end

;========================================
; conversion function
;========================================
to assignprop [partial-gene]
  let gl gene-length
  let SR1-coef (binary-to-decimal(graycode-to-binary(sublist partial-gene 0 (1 * gl)))) ;sensor and range
  ifelse Heterogeneous? [

    set UT-attraction (binary-to-decimal(graycode-to-binary(sublist partial-gene (1 * gl) (2 * gl))))
    set UP-attraction (binary-to-decimal(graycode-to-binary(sublist partial-gene (2 * gl) (3 * gl))))
    set UU-attraction (binary-to-decimal(graycode-to-binary(sublist partial-gene (3 * gl) (4 * gl))))
    set UU-repulsion (binary-to-decimal(graycode-to-binary(sublist partial-gene (4 * gl) (5 * gl))))
    set UU-attraction2 (binary-to-decimal(graycode-to-binary(sublist partial-gene (5 * gl) (6 * gl))))
    set UU-repulsion2 (binary-to-decimal(graycode-to-binary(sublist partial-gene (6 * gl) (7 * gl))))
    set UU-attraction3 (binary-to-decimal(graycode-to-binary(sublist partial-gene (7 * gl) (8 * gl))))
    set UU-repulsion3 (binary-to-decimal(graycode-to-binary(sublist partial-gene (8 * gl) (9 * gl))))
    set p-inc (binary-to-decimal(graycode-to-binary(sublist partial-gene (9 * gl) (10 * gl))))
    set p-coef (binary-to-decimal(graycode-to-binary(sublist partial-gene (10 * gl) (11 * gl)))) ; 40 44
    set p-coef p-coef / GBASE
  ]
  [

    set UT-attraction (binary-to-decimal(graycode-to-binary(sublist partial-gene (1 * gl) (2 * gl))))
    set UP-attraction (binary-to-decimal(graycode-to-binary(sublist partial-gene (2 * gl) (3 * gl))))
    set UU-attraction (binary-to-decimal(graycode-to-binary(sublist partial-gene (3 * gl) (4 * gl))))
    set UU-repulsion (binary-to-decimal(graycode-to-binary(sublist partial-gene (4 * gl) (5 * gl))))
    set p-inc (binary-to-decimal(graycode-to-binary(sublist partial-gene (5 * gl) (6 * gl))))
    set p-coef (binary-to-decimal(graycode-to-binary(sublist partial-gene (6 * gl) (7 * gl)))) ; 40 44
    set p-coef p-coef / GBASE
  ]
  ;set sensor-range (SR1-coef + 1) * MAX-sensor-range / SCM
  ;set uspeed (SCM - SR1-coef - 1) * MAX-uav-speed / SCM
  set sensor-range (SR1-coef ) * MAX-sensor-range / GBASE
  set uspeed (SCM - SR1-coef) * MAX-uav-speed / GBASE
  if sensor-range = 0 [set sensor-range 0.001]
end
to-report binary-to-decimal [bits]
  let bits-length (length bits) - 1
  let bits-decimal 0
  foreach bits [
    [i] -> set bits-decimal bits-decimal + ((2 ^ bits-length) * i)
    set bits-length bits-length - 1
  ]
  report bits-decimal
end

to-report graycode-to-binary [bits]
  let graycode bits
  let bits-index 0
  let binary n-values length bits [[]]

  foreach binary [
    ifelse bits-index = 0
    [set binary replace-item bits-index binary (item bits-index graycode)]
    [set binary replace-item bits-index binary (x-or (item bits-index graycode) (item (bits-index - 1) binary))]
    set bits-index bits-index + 1
  ]
    report binary
end

to-report x-or [a b]
  let c 0
  ifelse a = 1 xor b = 1
  [set c 1]
  [set c 0]
  report c
end

to-report index-of-min [index-list]
  report position (min index-list) index-list
end

;============================================
; PLOT
;============================================
to-report diversity
let distances []
let pop-list n-values population-size [[k] -> k]
foreach pop-list
[ [i]->
  let pop-list2 n-values (population-size - (i + 1)) [[j] -> j + i + 1]
  ;show pop-list2
foreach  pop-list2
[ [j]->
  set distances fput (hamming-distance (item j genetic-code-list) (item i genetic-code-list)) distances ] ]
let max-possible-distance-sum floor (population-size * population-size / 4)
report (sum distances) / max-possible-distance-sum
end

to-report hamming-distance [bits1 bits2]
  report (length remove true (map [ [b1 b2] -> b1 = b2 ] bits1 bits2)) / (length bits1)
end

to my-update-plots
let average-time ((sum generation-times) / (length generation-times))
set-current-plot "Fitness Plot"
set-current-plot-pen "Fitness" ;plot minimum-fitness
ifelse GA-mode = 1 [plot 1 - (min generation-times)][plot (max-run-time - (min generation-times)) / max-run-time ]

set-current-plot "Diversity Plot"
set-current-plot-pen "diversity"
plot diversity
show diversity
end
@#$#@#$#@
GRAPHICS-WINDOW
265
10
799
545
-1
-1
6.5
1
10
1
1
1
0
0
0
1
-40
40
-40
40
1
1
1
ticks
30.0

BUTTON
24
10
94
43
NIL
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SLIDER
24
49
196
82
population-size
population-size
0
100
12.0
1
1
NIL
HORIZONTAL

SLIDER
23
81
195
114
max-generation
max-generation
1
1000
1000.0
1
1
NIL
HORIZONTAL

BUTTON
98
10
161
43
NIL
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SLIDER
23
114
195
147
num-tests
num-tests
1
100
12.0
1
1
NIL
HORIZONTAL

SWITCH
846
409
1005
442
show-background?
show-background?
1
1
-1000

SLIDER
23
146
195
179
agent-population
agent-population
3
100
9.0
1
1
NIL
HORIZONTAL

SLIDER
23
179
195
212
target-population
target-population
1
100
65.0
1
1
NIL
HORIZONTAL

SLIDER
870
16
1070
49
MAX-sensor-range
MAX-sensor-range
1
100
9.0
1
1
patches
HORIZONTAL

SLIDER
23
276
195
309
uav-max-turn
uav-max-turn
0
100
20.0
1
1
NIL
HORIZONTAL

SLIDER
23
308
195
341
helis-speed
helis-speed
0
2
1.0
0.1
1
NIL
HORIZONTAL

SLIDER
871
113
1005
146
target-max-turn
target-max-turn
0
100
28.0
1
1
NIL
HORIZONTAL

SLIDER
871
80
1071
113
target-target-separation
target-target-separation
0
100
3.0
1
1
NIL
HORIZONTAL

SLIDER
1069
48
1241
81
target-sensor-range
target-sensor-range
0
100
11.0
1
1
NIL
HORIZONTAL

SLIDER
870
48
1071
81
target-uav-repulsion
target-uav-repulsion
0
100
17.0
1
1
NIL
HORIZONTAL

SLIDER
1069
80
1241
113
target-target-repulsion
target-target-repulsion
0
100
4.0
1
1
NIL
HORIZONTAL

SLIDER
23
341
195
374
MAX-uav-speed
MAX-uav-speed
0
2
1.0
0.1
1
NIL
HORIZONTAL

SLIDER
1070
16
1242
49
target-speed
target-speed
0
1
0.5
0.01
1
NIL
HORIZONTAL

SLIDER
23
211
195
244
max-run-time
max-run-time
0
10000
1146.0
1
1
NIL
HORIZONTAL

SLIDER
23
244
195
277
mutation-rate
mutation-rate
0
100
3.0
0.5
1
NIL
HORIZONTAL

SWITCH
1006
409
1134
442
show-energy
show-energy
1
1
-1000

BUTTON
870
146
1004
179
show best pop
let testt avg-time(best-pop-gene)\n;show best-pop-gene
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

INPUTBOX
1004
113
1241
173
save-file-name
hetero_speed_0dot5_gen.txt
1
0
String

PLOT
846
259
1074
409
Fitness Plot
generations
Fitness Value
0.0
20.0
0.0
1.0
true
false
"" ""
PENS
"Fitness" 1.0 0 -16777216 true "" ""

PLOT
1073
259
1278
409
Diversity Plot
Generations
Diversity
0.0
20.0
0.0
2.0
true
false
"" ""
PENS
"diversity" 1.0 0 -16777216 true "" ""

BUTTON
870
179
936
212
NIL
go-GA
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SWITCH
1131
409
1278
442
Heterogeneous?
Heterogeneous?
0
1
-1000

SWITCH
914
481
1071
514
Bayesian-method?
Bayesian-method?
0
1
-1000

SWITCH
1095
480
1233
513
test-bayesian?
test-bayesian?
1
1
-1000

SLIDER
23
374
195
407
MAX-pINC
MAX-pINC
0
100
16.0
1
1
NIL
HORIZONTAL

SLIDER
23
407
195
440
MAX-pDEC
MAX-pDEC
0
100
4.0
1
1
NIL
HORIZONTAL

BUTTON
1004
172
1156
205
save generation gene
file-open save-file-name\nfile-write genetic-code-list\nfile-close
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
1155
172
1273
205
load generation
file-open save-file-name\nset genetic-code-list file-read \nfile-close
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

@#$#@#$#@
## WHAT IS IT?

(a general understanding of what the model is trying to show or explain)

## HOW IT WORKS

(what rules the agents use to create the overall behavior of the model)

## HOW TO USE IT

(how to use the model, including a description of each of the items in the Interface tab)

## THINGS TO NOTICE

(suggested things for the user to notice while running the model)

## THINGS TO TRY

(suggested things for the user to try to do (move sliders, switches, etc.) with the model)

## EXTENDING THE MODEL

(suggested things to add or change in the Code tab to make the model more complicated, detailed, accurate, etc.)

## NETLOGO FEATURES

(interesting or unusual features of NetLogo that the model uses, particularly in the Code tab; or where workarounds were needed for missing features)

## RELATED MODELS

(models in the NetLogo Models Library and elsewhere which are of related interest)

## CREDITS AND REFERENCES

(a reference to the model's URL on the web if it has one, as well as any other necessary credits, citations, and links)
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.0
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@
