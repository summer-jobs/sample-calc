// A/B Testing Sample & Effect Size Calculator JS
// Fixed: Bonferroni correction now correctly impacts Binary N calculations.

// --- CONSTANTS ---
const POWER_ZB = 0.8416; // Z-score for 80% Power (1 - beta)
const VIF = 1.1; // Variance Inflation Factor
const EPSILON = 1e-9; 

// --- UI ELEMENT HELPERS ---
const $d = (id) => document.getElementById(id);

const inputs = {
 dataType: () => $d('dataTypeSelect').value,
 confidence: () => $d('confidenceSelect').value,
 groups: () => parseFloat($d('groups').value),
 baseline: () => {
  const rawValue = parseFloat($d('baseline').value);
  const dataType = inputs.dataType();
  if (isNaN(rawValue)) return NaN;
  if (dataType === 'binary') { return rawValue / 100; }
  return rawValue;
 },
 stdDev: () => parseFloat($d('stdDev').value),
 mde: () => {
  const inputVal = parseFloat($d('mde').value);
  const baseline = inputs.baseline();
  if (isNaN(inputVal) || isNaN(baseline)) return NaN;
  return baseline * (inputVal / 100); 
 },
 n: () => parseFloat($d('sampleSize').value)
};

// --- NORMAL INVERSE (Φ⁻¹) FUNCTION ---
function ppnd(p) {
    var a0 = 2.50662823884; var a1 = -18.61500062529; var a2 = 41.39119773534; var a3 = -25.44106049637;
    var b1 = -8.47351093090; var b2 = 23.08336743743; var b3 = -21.06224101826; var b4 = 3.13082909833;
    var c0 = -2.78718931138; var c1 = -2.29796479134; var c2 = 4.85014127135; var c3 = 2.32121276858;
    var d1 = 3.54388924762; var d2 = 1.63706781897; var r; var split = 0.42; var value;
    if ( Math.abs( p - 0.5 ) <= split ) {
        r = ( p - 0.5 ) * ( p - 0.5 );
        value = ( p - 0.5 ) * ( ( ( a3 * r + a2 ) * r + a1 ) * r + a0 ) / ( ( ( ( b4 * r + b3 ) * r + b2 ) * r + b1 ) * r + 1.0 );
    } else if ( 0.0 < p && p < 1.0 ) {
        if ( 0.5 < p ) { r = Math.sqrt ( - Math.log ( 1.0 - p ) ); } else { r = Math.sqrt ( - Math.log ( p ) ); }
        value = ( ( ( c3 * r + c2 ) * r + c1 ) * r + c0 ) / ( ( d2 * r + d1 ) * r + 1.0 );
        if ( p < 0.5 ) { value = - value; }
    } else { value = NaN; }
    return value;
}

// --- Z-alpha WITH MULTIPLE-GROUP (BONFERRONI) CORRECTION ---
function getZalpha() {
 const confStr = inputs.confidence(); const conf = parseFloat(confStr);
 if (isNaN(conf) || conf <= 0 || conf >= 1) { return 1.96; }
 // Comparisons = Groups - 1
 const mComparisons = (!isNaN(inputs.groups()) && inputs.groups() > 2) ? (inputs.groups() - 1) : 1;
 if (mComparisons === 1) { 
  const fixedScores = { '0.90': 1.645, '0.95': 1.960, '0.99': 2.576 }; 
  return fixedScores[confStr] || 1.960; 
 }
 const alphaPerComp = (1 - conf) / mComparisons;
 return ppnd(1 - (alphaPerComp / 2));
}

// --- UPDATED BINARY N FORMULA (Accepts Zalpha directly) ---
function calculateBinaryN(Zalpha, power_level, p_baseline, delta_absolute_MDE) {
    if (p_baseline > 0.5) {
        p_baseline = 1.0 - p_baseline;
        delta_absolute_MDE = -delta_absolute_MDE;
    }
    
    var t_alpha2 = Zalpha; // Now using the corrected Z-score passed from getZalpha()
    var t_beta = ppnd(power_level);

    var sd1 = Math.sqrt(2 * p_baseline * (1.0 - p_baseline));
    var sd2 = Math.sqrt(p_baseline * (1.0 - p_baseline) + (p_baseline + delta_absolute_MDE) * (1.0 - p_baseline - delta_absolute_MDE));

    return (t_alpha2 * sd1 + t_beta * sd2) * (t_alpha2 * sd1 + t_beta * sd2) / (delta_absolute_MDE * delta_absolute_MDE);
}

// --- STATISTICAL CALCULATION FUNCTIONS ---

function calculateN(Zalpha, MDE) {
 const P_A = inputs.baseline();
 const dataType = inputs.dataType();
 const Z_sum_sq = Math.pow(Zalpha + POWER_ZB, 2);
 let N_raw;

 if (dataType === 'binary') {
  // Pass the corrected Zalpha directly into the binary formula
  N_raw = calculateBinaryN(Zalpha, 0.80, P_A, MDE); 
 } else {
  const sigma = inputs.stdDev();
  N_raw = (2 * Math.pow(sigma, 2) * Z_sum_sq) / Math.pow(MDE, 2);
 }
 return N_raw * VIF;
}

function calculateMDE(Zalpha, N_target) {
 const P_A = inputs.baseline();
 const dataType = inputs.dataType();
 const N_target_raw = N_target / VIF;
 const Z_sum = Zalpha + POWER_ZB;
 
 if (dataType === 'binary') {
  const pooled_variance = (P_A * (1 - P_A) + (P_A * (1 - P_A))) / 2;
  return Z_sum * Math.sqrt((2 * pooled_variance) / N_target_raw);
 } else {
  const sigma = inputs.stdDev();
  return Math.sqrt((2 * Math.pow(sigma, 2) * Math.pow(Z_sum, 2)) / N_target_raw);
 }
}

// --- UI UPDATES & HANDLERS ---
let isCalculating = false;
let lastInputSource = 'mde'; 

function handleInput(source) {
 if (isCalculating) return;
 if (source === 'n') { $d('mde').value = ""; } else { $d('sampleSize').value = ""; }
 lastInputSource = source;
 calculate();
}

function updateUI() {
 const isContinuous = inputs.dataType() === 'continuous';
 if (isContinuous) { $d('sdContainer').classList.remove('hidden'); } else { $d('sdContainer').classList.add('hidden'); }
 $d('baselineLabel').textContent = isContinuous ? 'Baseline Mean' : 'Baseline Rate (%)';
 $d('mdeLabel').textContent = 'Minimum Detectable Relative Lift (%)';
 $d('mde').value = ""; $d('sampleSize').value = ""; $d('statusMessageOutput').textContent = '';
 $d('baseline').value = ""; $d('stdDev').value = "";
 $d('totalUsersNeeded').textContent = '—';
 lastInputSource = 'mde'; 
 calculate();
}

function calculate() {
 isCalculating = true;
 const Zalpha = getZalpha(); 
 const groups = inputs.groups();
 const P_A = inputs.baseline();
 const sigma = inputs.stdDev();
 let nPerGroup, mdeAbsolute;
 let statusMessage = '';
 
 const isContinuous = inputs.dataType() === 'continuous';
 if (isContinuous && (isNaN(sigma) || sigma <= 0)) {
  statusMessage = 'ERROR: Please enter a valid Standard Deviation (> 0) for continuous data.';
 } else {
  const nInput = inputs.n();
  const mdeInput = inputs.mde();
  if (lastInputSource === 'n') {
   if (nInput > 0 && P_A > 0) {
    mdeAbsolute = calculateMDE(Zalpha, nInput);
    nPerGroup = nInput;
   }
  } else {
   if (mdeInput > 0 && P_A > 0) {
    nPerGroup = calculateN(Zalpha, mdeInput);
    mdeAbsolute = mdeInput;
   }
  }
 }

 const sampleSizeInput = $d('sampleSize');
 const mdeInputUI = $d('mde');
 const totalUsersNeededOutput = $d('totalUsersNeeded');
 
 if (statusMessage === '' && !isNaN(nPerGroup) && !isNaN(mdeAbsolute)) {
  if (lastInputSource !== 'n') {
   sampleSizeInput.value = Math.ceil(nPerGroup);
  } else {
   const mdeRelativePercent = (mdeAbsolute / P_A) * 100;
   mdeInputUI.value = mdeRelativePercent.toFixed(2);
  }
  totalUsersNeededOutput.textContent = Math.ceil(Math.ceil(nPerGroup) * groups);
 }
 isCalculating = false;
}
