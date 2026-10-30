package io.github.motohiroogawa.msentity.similarity;

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

  /** cosine, 0.01 Da tolerance, ppm disabled, exponent 1, all threads. */
  public static Options defaults() {
    return new Options(Method.COSINE, 0.01, 0.0, 1.0, 0);
  }

  /**
   * Compatibility constructor: the former bin width is used as the Da tolerance.
   *
   * @deprecated use the constructor with {@code mzToleranceDa} and {@code mzTolerancePpm}
   */
  @Deprecated
  public Options(Method method, double binWidth, double intensityExponent, int threads) {
    this(method, binWidth, 0.0, intensityExponent, threads);
  }

  public Options withMethod(Method value) {
    return new Options(value, mzToleranceDa, mzTolerancePpm, intensityExponent, threads);
  }

  /**
   * Former name of {@link #withMzToleranceDa(double)}.
   *
   * @deprecated use {@link #withMzToleranceDa(double)}
   */
  @Deprecated
  public Options withBinWidth(double value) {
    return withMzToleranceDa(value);
  }

  public Options withMzToleranceDa(double value) {
    return new Options(method, value, mzTolerancePpm, intensityExponent, threads);
  }

  public Options withMzTolerancePpm(double value) {
    return new Options(method, mzToleranceDa, value, intensityExponent, threads);
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
