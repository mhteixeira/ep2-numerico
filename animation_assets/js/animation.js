let scl = 6;
let amplitude = 3000;
let mode = 0  ;

const initialTime = 0;
const endTime = 6;
const dt = (3 * (endTime - initialTime)) ** -1;
// let dt = 1/4
let t = initialTime;

let joints = [];
let nosTotais;
let posicoes;
let matrizDeBarras;
let x = [];

function setup() {

  createCanvas(400, 400);
  frameRate(30);
  
  sliderAmplitude = createSlider(50, 250, 50, 10);
  sliderAmplitude.position(270, 368);
  sliderAmplitude.style('width', '80px');

  sliderMode = createSlider(0, 4, 0, 1);
  sliderMode.position(160, 25);
  sliderMode.style('width', '80px');

  infos = loadTrussesData(inputC);
  nosTotais = infos[0];
  posicao_x = infos[3];
  posicao_y = infos[4];
  matrizDeBarras = infos[5];


  for (i = 0; i <= nosTotais; i++) {
    joints[i] = new Joint(scl * posicao_x[i], scl * posicao_y[i], i);
  }

  createP(`Animação das treliças vibrando para o EP2 <br/>
            Para mudar tanto o modo de vibração em análise quanto a amplitude do movimento, use os sliders`);

}

function draw() {
  background(51);
  let mode = sliderMode.value();
  frequency = getFrequency(splitLines(frequencias), mode);
  vibrationMode = getVibrationMode(splitLines(modosDeVibracao), mode);
  let amplitude = sliderAmplitude.value();
  translate(80, height-30);
  scale(1, -1);
  
  
  
  //Desenho do chão
  fill(110);
  rect(0, 0, -700, 80);

  //Desenho das barras
  for (j = 0; j < barras; j++) {
    [primeiroNo, segundoNo] = matrizDeBarras[j];
    stroke(255);
    strokeWeight(2);
    line(
      joints[primeiroNo].x,
      joints[primeiroNo].y,
      joints[segundoNo].x,
      joints[segundoNo].y
    );
  }

  // Desenho dos nós
  for (k = 1; k <= nosTotais; k++)
    joints[k].show();

  angleMode(RADIANS)

  x = calculateX(vibrationMode, frequency, t);
  // Deslocamentos
  for (k = 1; k <= nosTotais - 2; k++) {
      joints[k].update(
        scl * amplitude * x[ 2*(k - 1) ], 
        scl * amplitude * x[ 2*(k - 1) + 1 ]
      );
  }

  fill(255);
  textSize(12);
  textStyle(NORMAL);
  textAlign(LEFT, CENTER);
  push();
  translate(width/2, height/2)
  rotate(PI);
  scale(-1, 1);
  text("Modo de vibração: " + (mode + 1), -260, -140);
  text("Tempo: " + t.toFixed(3) + "s", -260, -110);
  text("Deslocamentos ampliados em " + amplitude + " vezes", -260, 200);
  pop();

  t += dt;
  if (t >= endTime) t = initialTime;
}
