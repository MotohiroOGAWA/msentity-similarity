package io.github.motohiroogawa.msentity.similarity;

import static io.github.motohiroogawa.msentity.similarity.Native.call;
import static io.github.motohiroogawa.msentity.similarity.Native.check;

import java.lang.foreign.Arena;
import java.lang.foreign.MemorySegment;
import java.lang.foreign.ValueLayout;

/**
 * Score settings. Immutable; change values with the {@code with…} methods:
 *
 * <pre>{@code
 * Options options = Options.defaults().withMethod(Method.REVERSE_COSINE).withMzTolerancePpm(10);
 * }</pre>
 *
 * @param method cosine or reverse cosine
 * @param mzToleranceDa absolute m/z tolerance in Da; 0 disables it
 * @param mzTolerancePpm symmetric m/z tolerance in ppm; 0 disables it
 * @param intensityExponent intensities are raised to this power (&gt; 0) before scoring
 * @param threads worker threads; 0 = all hardware threads
 */
public record Options(Method method, double mzToleranceDa, double mzTolerancePpm,
                      double intensityExponent, int threads) {
  /** The scoring method. */
  public enum Method {
    COSINE,
    REVERSE_COSINE
  }

  /** cosine, 0.01 Da or 10 ppm tolerance, exponent 1, all threads. */
  public static Options defaults() {
    return new Options(Method.COSINE, 0.01, 10.0, 1.0, 0);
  }

  public Options withMethod(Method value) {
    return new Options(value, mzToleranceDa, mzTolerancePpm, intensityExponent, threads);
  }

  public Options withMzToleranceDa(double value) {
    return new Options(method, value, mzTolerancePpm, intensityExponent, threads);
  }

  public Options withMzTolerancePpm(double value) {
    return new Options(method, mzToleranceDa, value, intensityExponent, threads);
  }

  /**
   * Both tolerances from text such as {@code "0.01Da,10ppm"} or {@code "10ppm"}: {@code <value>Da}
   * and/or {@code <value>ppm} in either order, units case-insensitive; a unit left out becomes 0.
   */
  public Options withMzTolerance(String text) {
    try (Arena arena = Arena.ofConfined()) {
      MemorySegment options = toNative(arena);
      MemorySegment t = Native.string(arena, text);
      check(call(() -> (int) Native.parseMzTolerance.invoke(t, options)));
      return new Options(method, options.get(ValueLayout.JAVA_DOUBLE, 8),
                         options.get(ValueLayout.JAVA_DOUBLE, 16), intensityExponent, threads);
    }
  }

  /** Both tolerances as canonical text, e.g. {@code "0.01Da,10ppm"} (a zero term is left out). */
  public String formatMzTolerance() {
    try (Arena arena = Arena.ofConfined()) {
      MemorySegment options = toNative(arena);
      MemorySegment slot = arena.allocate(ValueLayout.ADDRESS);
      check(call(() -> (int) Native.formatMzTolerance.invoke(options, slot)));
      return Native.takeString(slot);
    }
  }

  public Options withIntensityExponent(double value) {
    return new Options(method, mzToleranceDa, mzTolerancePpm, value, threads);
  }

  public Options withThreads(int value) {
    return new Options(method, mzToleranceDa, mzTolerancePpm, intensityExponent, value);
  }

  /** struct mssim_options in {@code arena}. */
  MemorySegment toNative(Arena arena) {
    MemorySegment out = arena.allocate(Native.OPTIONS);
    out.set(ValueLayout.JAVA_INT, 0, method == Method.REVERSE_COSINE ? 1 : 0);
    out.set(ValueLayout.JAVA_INT, 4, threads);
    out.set(ValueLayout.JAVA_DOUBLE, 8, mzToleranceDa);
    out.set(ValueLayout.JAVA_DOUBLE, 16, mzTolerancePpm);
    out.set(ValueLayout.JAVA_DOUBLE, 24, intensityExponent);
    return out;
  }
}
