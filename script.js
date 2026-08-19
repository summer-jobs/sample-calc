// A/B Testing Sample & Effect Size Calculator
//
// Method: Fleiss-style normal approximation for binary outcomes (null-hypothesis
// SD pooled at the average of the two arms, alternative-hypothesis SD unpooled),
// and a normal approximation with a small-sample df adjustment for continuous
// outcomes. Bonferroni correction on alpha when there are more than two groups.
// A variance inflation factor is applied last, for conservative planning.

// --- CONSTANTS ---
const POWER = 0.80;          // Fixed statistical power (1 - beta)
const VIF = 1.1;             // Variance inflation factor, applied to the final N
const T_DF_ADJUSTMENT = 2;   // Continuous only: approximates t vs. z (Guenther/Snedecor)
const MDE_SOLVER_ITERATIONS = 200;

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
// Beasley-Springer (AS 111). Accurate to ~1e-9 across the range used here.
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

// Single source of truth for the power z-score.
const Z_BETA = ppnd(POWER);

// --- Z-alpha WITH MULTIPLE-GROUP (BONFERRONI) CORRECTION ---
// Comparisons are counted against a single control arm, so m = groups - 1.
function getZalpha() {
 const conf = parseFloat(inputs.confidence());
 if (isNaN(conf) || conf <= 0 || conf >= 1) { return ppnd(0.975); }
 const groups = inputs.groups();
 const mComparisons = (!isNaN(groups) && groups > 2) ? (groups - 1) : 1;
 const alphaPerComparison = (1 - conf) / mComparisons;
 return ppnd(1 - (alphaPerComparison / 2));
}

// --- BINARY N ---
// N per group, before VIF. p2 = p1 + delta.
// The null-hypothesis SD (sd1) uses the AVERAGE of the two arms, because under
// the null both arms sit at that common rate. Using p1 here understates the
// required N, and the error grows with the size of the lift.
function calculateBinaryN(Zalpha, Zbeta, p1, delta) {
    const p2 = p1 + delta;
    const pBar = p1 + delta / 2;

    const sd1 = Math.sqrt(2 * pBar * (1.0 - pBar));
    const sd2 = Math.sqrt(p1 * (1.0 - p1) + p2 * (1.0 - p2));

    const numerator = Zalpha * sd1 + Zbeta * sd2;
    return (numerator * numerator) / (delta * delta);
}

// --- CONTINUOUS N ---
function calculateContinuousN(Zalpha, Zbeta, sigma, delta) {
    const zSum = Zalpha + Zbeta;
    return (2 * sigma * sigma * zSum * zSum) / (delta * delta);
}

// --- FORWARD: MDE -> N per group ---
function calculateN(Zalpha, MDE) {
 if (inputs.dataType() === 'binary') {
  return calculateBinaryN(Zalpha, Z_BETA, inputs.baseline(), MDE) * VIF;
 }
 return (calculateContinuousN(Zalpha, Z_BETA, inputs.stdDev(), MDE) + T_DF_ADJUSTMENT) * VIF;
}

// --- REVERSE: N per group -> MDE (absolute) ---
// Binary has no clean closed form, so it is solved numerically against the same
// function used in the forward direction. This keeps the two fields consistent:
// a value typed into either box round-trips back to the other.
function calculateMDE(Zalpha, N_target) {
 const dataType = inputs.dataType();

 if (dataType !== 'binary') {
  const rawTarget = (N_target / VIF) - T_DF_ADJUSTMENT;
  if (rawTarget <= 0) { return NaN; }
  const sigma = inputs.stdDev();
  const zSum = Zalpha + Z_BETA;
  return Math.sqrt((2 * sigma * sigma * zSum * zSum) / rawTarget);
 }

 const p1 = inputs.baseline();
 let low = 0;
 let high = (1 - p1) * 0.999999; // keep p2 below 100%

 // N falls as the lift grows. If even the largest possible lift needs more
 // users than are available, there is no solution to find.
 if (calculateBinaryN(Zalpha, Z_BETA, p1, high) * VIF > N_target) { return NaN; }

 for (let i = 0; i < MDE_SOLVER_ITERATIONS; i++) {
  const mid = (low + high) / 2;
  if (mid <= 0) { break; }
  if (calculateBinaryN(Zalpha, Z_BETA, p1, mid) * VIF > N_target) {
   low = mid;
  } else {
   high = mid;
  }
 }
 return high;
}

// --- VALIDATION ---
// Returns an empty string when the inputs are usable, or a message naming the
// field to fix.
function validateInputs(source) {
 const dataType = inputs.dataType();
 const groups = inputs.groups();
 const baseline = inputs.baseline();

 if (isNaN(groups) || groups < 2 || groups !== Math.floor(groups)) {
  return 'Enter a whole number of groups (2 or more).';
 }

 if (dataType === 'binary') {
  if (isNaN(baseline)) { return ''; }
  if (baseline <= 0 || baseline >= 1) {
   return 'Baseline rate must be between 0 and 100 (exclusive).';
  }
 } else {
  if (!isNaN(baseline) && baseline <= 0) {
   return 'Baseline mean must be greater than 0, since the lift is relative to it.';
  }
  const sigma = inputs.stdDev();
  if (isNaN(sigma) && !isNaN(baseline)) { return 'Enter a standard deviation for continuous data.'; }
  if (!isNaN(sigma) && sigma <= 0) { return 'Standard deviation must be greater than 0.'; }
 }

 if (source === 'mde') {
  const mdeRelative = parseFloat($d('mde').value);
  if (!isNaN(mdeRelative) && mdeRelative <= 0) {
   return 'Minimum detectable lift must be greater than 0.';
  }
  if (dataType === 'binary' && !isNaN(baseline) && !isNaN(mdeRelative)) {
   if (baseline * (1 + mdeRelative / 100) >= 1) {
    const maxLift = ((1 - baseline) / baseline) * 100;
    return 'That lift puts the treatment rate at or above 100%. With a ' +
           (baseline * 100).toFixed(1) + '% baseline, the largest possible lift is ' +
           maxLift.toFixed(1) + '%.';
   }
  }
 } else {
  const n = inputs.n();
  if (!isNaN(n) && n < 2) { return 'Sample size per group must be at least 2.'; }
 }

 return '';
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
 $d('mde').value = ""; $d('sampleSize').value = "";
 $d('baseline').value = ""; $d('stdDev').value = "";
 $d('statusMessageOutput').textContent = '';
 $d('totalUsersNeeded').textContent = '—';
 lastInputSource = 'mde';
 calculate();
}

function calculate() {
 isCalculating = true;

 const statusOutput = $d('statusMessageOutput');
 const totalOutput = $d('totalUsersNeeded');
 const sampleSizeInput = $d('sampleSize');
 const mdeInputUI = $d('mde');

 const statusMessage = validateInputs(lastInputSource);
 if (statusMessage !== '') {
  statusOutput.textContent = statusMessage;
  totalOutput.textContent = '—';
  if (lastInputSource === 'mde') { sampleSizeInput.value = ""; } else { mdeInputUI.value = ""; }
  isCalculating = false;
  return;
 }
 statusOutput.textContent = '';

 const Zalpha = getZalpha();
 const groups = inputs.groups();
 const baseline = inputs.baseline();
 let nPerGroup = NaN;
 let mdeAbsolute = NaN;

 if (lastInputSource === 'n') {
  const nInput = inputs.n();
  if (nInput > 0 && baseline > 0) {
   mdeAbsolute = calculateMDE(Zalpha, nInput);
   nPerGroup = nInput;
  }
 } else {
  const mdeInput = inputs.mde();
  if (mdeInput > 0 && baseline > 0) {
   nPerGroup = calculateN(Zalpha, mdeInput);
   mdeAbsolute = mdeInput;
  }
 }

 if (!isFinite(nPerGroup) || !isFinite(mdeAbsolute)) {
  if (lastInputSource === 'n') {
   mdeInputUI.value = "";
   if (inputs.n() > 0 && baseline > 0) {
    statusOutput.textContent = 'That sample size is too small to detect any lift at this baseline and confidence level.';
   }
  } else {
   sampleSizeInput.value = "";
  }
  totalOutput.textContent = '—';
  isCalculating = false;
  return;
 }

 const nPerGroupRounded = Math.ceil(nPerGroup);
 if (lastInputSource === 'mde') {
  sampleSizeInput.value = nPerGroupRounded;
 } else {
  mdeInputUI.value = ((mdeAbsolute / baseline) * 100).toFixed(2);
 }
 totalOutput.textContent = (nPerGroupRounded * groups).toLocaleString();

 isCalculating = false;
}
