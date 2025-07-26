const width = 100;
const height = 100;
const dt = 1 / 60;
let diffusionRate = 0.2;
let diffusionIterations = 10;
let initialVelocity = 10;

class GridSpace {
  constructor(velocityX = 0, velocityY = 0, density = 0) {
    this.velocityX = velocityX;
    this.velocityY = velocityY;
    this.density = density;
  }
}

const grid = Array.from({ length: height }, () =>
  Array.from({ length: width }, () => new GridSpace())
);

function lerp(a, b, t) {
  return a * (1 - t) + b * t;
}

function diffuse(diffusionRate, iterations) {
  const k = diffusionRate * dt;
  const newGrid = Array(height)
    .fill(null)
    .map(() => Array(width).fill(0));

  for (let y = 0; y < height; y++)
    for (let x = 0; x < width; x++) newGrid[y][x] = grid[y][x].density;

  //needed to solve the system of equations
  for (let iter = 0; iter < iterations; iter++) {
    for (let y = 0; y < height; y++) {
      for (let x = 0; x < width; x++) {
        const left = x === 0 ? newGrid[y][x] : newGrid[y][x - 1];
        const right = x === width - 1 ? newGrid[y][x] : newGrid[y][x + 1];
        const up = y === 0 ? newGrid[y][x] : newGrid[y - 1][x];
        const down = y === height - 1 ? newGrid[y][x] : newGrid[y + 1][x];

        const neighborsAvg = (left + right + up + down) / 4;

        newGrid[y][x] = (newGrid[y][x] + k * neighborsAvg) / (1 + k);
      }
    }
  }

  for (let y = 0; y < height; y++)
    for (let x = 0; x < width; x++) grid[y][x].density = newGrid[y][x];
}

//get rid of the divergence of the vector field
function projectVelocity(grid, width, height, iterations = 40) {
  const div = Array.from({ length: height }, () => Array(width).fill(0));
  const phi = Array.from({ length: height }, () => Array(width).fill(0));

  for (let y = 1; y < height - 1; y++) {
    for (let x = 1; x < width - 1; x++) {
      div[y][x] =
        0.5 *
        (grid[y][x + 1].velocityX -
          grid[y][x - 1].velocityX +
          grid[y + 1][x].velocityY -
          grid[y - 1][x].velocityY);
      phi[y][x] = 0;
    }
  }

  for (let iter = 0; iter < iterations; iter++) {
    for (let y = 1; y < height - 1; y++) {
      for (let x = 1; x < width - 1; x++) {
        phi[y][x] =
          (phi[y][x + 1] +
            phi[y][x - 1] +
            phi[y + 1][x] +
            phi[y - 1][x] -
            div[y][x]) /
          4;
      }
    }
  }

  for (let y = 1; y < height - 1; y++) {
    for (let x = 1; x < width - 1; x++) {
      grid[y][x].velocityX -= 0.5 * (phi[y][x + 1] - phi[y][x - 1]);
      grid[y][x].velocityY -= 0.5 * (phi[y + 1][x] - phi[y - 1][x]);
    }
  }

  for (let x = 0; x < width; x++) {
    grid[0][x].velocityX = 0;
    grid[0][x].velocityY = 0;
    grid[height - 1][x].velocityX = 0;
    grid[height - 1][x].velocityY = 0;
  }
  for (let y = 0; y < height; y++) {
    grid[y][0].velocityX = 0;
    grid[y][0].velocityY = 0;
    grid[y][width - 1].velocityX = 0;
    grid[y][width - 1].velocityY = 0;
  }
}

function advectDensity() {
  const newGrid = Array(height)
    .fill(null)
    .map(() => Array(width).fill(0));

  for (let y = 0; y < height; y++)
    for (let x = 0; x < width; x++) newGrid[y][x] = grid[y][x].density;

  for (let y = 0; y < height; y++) {
    for (let x = 0; x < width; x++) {
      let fx = x - grid[y][x].velocityX * dt;
      let fy = y - grid[y][x].velocityY * dt;

      fx = Math.min(Math.max(fx, 0), width - 1);
      fy = Math.min(Math.max(fy, 0), height - 1);

      let ix = Math.floor(fx);
      let iy = Math.floor(fy);
      let jx = fx - ix;
      let jy = fy - iy;

      const getDensity = (yy, xx) =>
        grid[Math.min(height - 1, Math.max(0, yy))][
          Math.min(width - 1, Math.max(0, xx))
        ].density;

      const valTopLeft = getDensity(iy, ix);
      const valTopRight = getDensity(iy, ix + 1);
      const valBottomLeft = getDensity(iy + 1, ix);
      const valBottomRight = getDensity(iy + 1, ix + 1);

      const valTop = lerp(valTopLeft, valTopRight, jx);
      const valBot = lerp(valBottomLeft, valBottomRight, jx);

      newGrid[y][x] = lerp(valTop, valBot, jy);
    }
  }

  for (let y = 0; y < height; y++)
    for (let x = 0; x < width; x++) grid[y][x].density = newGrid[y][x];
}

function advectVelocity() {
  const newGrid = Array(height)
    .fill(null)
    .map(() =>
      Array(width)
        .fill(null)
        .map(() => ({ velocityX: 0, velocityY: 0 }))
    );

  for (let y = 0; y < height; y++) {
    for (let x = 0; x < width; x++) {
      let fx = x - grid[y][x].velocityX * dt;
      let fy = y - grid[y][x].velocityY * dt;

      let ix = Math.floor(fx);
      let iy = Math.floor(fy);
      let jx = fx - ix;
      let jy = fy - iy;

      const getComp = (yy, xx, comp) =>
        grid[Math.min(height - 1, Math.max(0, yy))][
          Math.min(width - 1, Math.max(0, xx))
        ][comp];

      const vTopLeftX = getComp(iy, ix, "velocityX");
      const vTopRightX = getComp(iy, ix + 1, "velocityX");
      const vBottomLeftX = getComp(iy + 1, ix, "velocityX");
      const vBottomRightX = getComp(iy + 1, ix + 1, "velocityX");

      const vTopX = lerp(vTopLeftX, vTopRightX, jx);
      const vBottomX = lerp(vBottomLeftX, vBottomRightX, jx);
      const interpolatedX = lerp(vTopX, vBottomX, jy);

      const vTopLeftY = getComp(iy, ix, "velocityY");
      const vTopRightY = getComp(iy, ix + 1, "velocityY");
      const vBottomLeftY = getComp(iy + 1, ix, "velocityY");
      const vBottomRightY = getComp(iy + 1, ix + 1, "velocityY");

      const vTopY = lerp(vTopLeftY, vTopRightY, jx);
      const vBottomY = lerp(vBottomLeftY, vBottomRightY, jx);
      const interpolatedY = lerp(vTopY, vBottomY, jy);

      newGrid[y][x].velocityX = interpolatedX;
      newGrid[y][x].velocityY = interpolatedY;
    }
  }

  for (let y = 0; y < height; y++)
    for (let x = 0; x < width; x++) {
      grid[y][x].velocityX = newGrid[y][x].velocityX;
      grid[y][x].velocityY = newGrid[y][x].velocityY;
    }
}

const canvas = document.getElementById("fluidCanvas");
const ctx = canvas.getContext("2d");
const imageData = ctx.createImageData(width, height);
const pixels = imageData.data;

function render() {
  for (let y = 0; y < height; y++) {
    for (let x = 0; x < width; x++) {
      const i = (y * width + x) * 4;
      const density = Math.min(1, grid[y][x].density);
      const value = Math.floor(density * 255);
      pixels[i] = value;
      pixels[i + 1] = value;
      pixels[i + 2] = value;
      pixels[i + 3] = 255;
    }
  }
  ctx.putImageData(imageData, 0, 0);
}

let isMouseDown = false;

const velocitySlider = document.getElementById("velocitySlider");
const velocityValueSpan = document.getElementById("velocityValue");

if (velocitySlider && velocityValueSpan) {
  velocityValueSpan.textContent = velocitySlider.value;

  velocitySlider.addEventListener("input", (e) => {
    initialVelocity = parseFloat(e.target.value);
    velocityValueSpan.textContent = initialVelocity;
  });
}

const iterationSlider = document.getElementById("iterationSlider");
const iterationSliderSpan = document.getElementById("iterationValue");

if (iterationSlider && iterationSliderSpan) {
  iterationSliderSpan.textContent = iterationSlider.value;

  iterationSlider.addEventListener("input", (e) => {
    diffusionIterations = parseFloat(e.target.value);
    iterationSliderSpan.textContent = diffusionIterations;
  });
}

const diffusionSlider = document.getElementById("diffusionSlider");
const diffusionValueSpan = document.getElementById("diffusionValue");

if (diffusionSlider && diffusionValueSpan) {
  diffusionValueSpan.textContent = diffusionSlider.value;

  diffusionSlider.addEventListener("input", (e) => {
    diffusionRate = parseFloat(e.target.value);
    diffusionValueSpan.textContent = diffusionRate;
  });
}

canvas.addEventListener("mousedown", () => (isMouseDown = true));
canvas.addEventListener("mouseup", () => (isMouseDown = false));
canvas.addEventListener("mouseleave", () => (isMouseDown = false));
canvas.addEventListener("mousemove", (e) => {
  if (!isMouseDown) return;

  const rect = canvas.getBoundingClientRect();
  const x = Math.floor(((e.clientX - rect.left) * width) / rect.width);
  const y = Math.floor(((e.clientY - rect.top) * height) / rect.height);

  const radius = 2;
  const chance = 0.8;

  for (let dy = -radius; dy <= radius; dy++) {
    for (let dx = -radius; dx <= radius; dx++) {
      const nx = x + dx;
      const ny = y + dy;

      if (
        nx > 0 &&
        nx < width - 1 &&
        ny > 0 &&
        ny < height - 1 &&
        dx * dx + dy * dy <= radius * radius
      ) {
        if (Math.random() < chance) {
          grid[ny][nx].density = 1;

          if (
            typeof canvas.prevX === "number" &&
            typeof canvas.prevY === "number"
          ) {
            let vx = nx - canvas.prevX + (Math.random() - 0.5);
            let vy = ny - canvas.prevY + (Math.random() - 0.5);

            vx = Math.max(-1, Math.min(2, vx));
            vy = Math.max(-1, Math.min(2, vy));

            grid[ny][nx].velocityX += vx * initialVelocity;
            grid[ny][nx].velocityY += vy * initialVelocity;
          }
        }
      }
    }
  }

  canvas.prevX = x;
  canvas.prevY = y;
});

const clearBtn = document.getElementById("clearCanvas");
clearBtn.addEventListener("click", () => {
  for (let y = 0; y < height; y++)
    for (let x = 0; x < width; x++) grid[y][x] = new GridSpace();
});

function updateFluid() {
  projectVelocity(grid, width, height);

  advectVelocity();

  projectVelocity(grid, width, height);

  diffuse(diffusionRate, diffusionIterations);

  advectDensity();
}

// Main loop
function loop() {
  updateFluid();
  render();
  requestAnimationFrame(loop);
}

loop();
